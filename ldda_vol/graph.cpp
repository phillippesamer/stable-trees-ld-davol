#include "graph.h"

Graph::Graph()
{
    using_matrix = false;
    using_lemon = false;
    lemon_graph_modified = false;
    num_vertices = 0;
    num_edges = 0;
    mst_weight = 0;
}

Graph::Graph(long n, long m)
{
    using_matrix = false;
    using_lemon = false;
    lemon_graph_modified = false;
    num_vertices = n;
    num_edges = m;
    mst_weight = numeric_limits<double>::max();   // flag for: mst not computed

    adj_list.reserve(n);
    adj_list.insert( adj_list.begin(), n, list<long>() );

    s.reserve(m);
    t.reserve(m);
    w.reserve(m);
    mst_vector.reserve(m);
}

Graph::~Graph()
{
    adj_list.clear();

    s.clear();
    t.clear();
    w.clear();
    mst_vector.clear();

	if (using_matrix)
        free_index_matrix();

    if (using_lemon)
    {
        lemon_vertices.clear();
        lemon_edges.clear();
        delete lemon_weight;
        delete lemon_edges_inverted_index;
        delete lemon_graph;
        delete opposite_weights;
    }
}

void Graph::init_index_matrix()
{
	using_matrix = true;

    index_matrix = new long*[num_vertices];
    for (long i=0; i<num_vertices; ++i)
    {
        index_matrix[i] = new long[num_vertices];
        for (long j=0; j<num_vertices; ++j)
            index_matrix[i][j] = -1;
    }
}

void Graph::free_index_matrix()
{
    for (long i = 0; i<num_vertices; ++i)
        delete[] index_matrix[i];

    delete[] index_matrix;
}

void Graph::init_lemon()
{
    using_lemon = true;

    lemon_graph = new ListGraph();

    lemon_vertices.reserve(num_vertices);
    for (long i=0; i<num_vertices; ++i)
        lemon_vertices.push_back(lemon_graph->addNode());

    lemon_edges.reserve(num_edges);
    lemon_weight = new ListGraph::EdgeMap<double>(*lemon_graph);
    lemon_edges_inverted_index = new ListGraph::EdgeMap<long>(*lemon_graph);

    // volume only
    opposite_weights = new ListGraph::EdgeMap<double>(*lemon_graph);
}

void Graph::update_single_weight(long idx, double new_weight)
{
    // weight in the edge list
    this->w[idx] = new_weight;

    // weight in lemon's adjacency list
    ListGraph::Edge e = this->lemon_edges[idx];
    (*lemon_weight)[e] = new_weight;
}

void Graph::update_all_weights(vector<double> new_weights)
{
    // weight in the edge list
    this->w.clear();
    this->w = vector<double>(new_weights);

    // weight in lemon's adjacency list
    for (long i=0; i<num_edges; ++i)
    {
        ListGraph::Edge e = this->lemon_edges[i];
        (*lemon_weight)[e] = new_weights[i];
    }
}

void Graph::lemon_delete_edge(long delete_index)
{
    this->lemon_graph_modified = true;

    this->lemon_graph->erase(this->lemon_edges[delete_index]);
}

vector<long> Graph::lemon_contract_edge(long contract_index)
{
    /// Contracts edge and returns the indices of any parallel edge dropped

    this->lemon_graph_modified = true;

    long vertex_1 = this->s[contract_index];
    long vertex_2 = this->t[contract_index];

    // return value: indices IN THE ORIGINAL GRAPH of parallel edges dropped
    vector<long> dropped_indices = vector<long>();

    vector<ListGraph::Edge> parallel_edges = 
        lemon_parallel_edges_if_contract( *lemon_graph,
                                          *lemon_weight,
                                           lemon_vertices[vertex_1],
                                           lemon_vertices[vertex_2] );

    // references to make notation less cumbersome
    ListGraph &g = *lemon_graph;
    ListGraph::Node &x = lemon_vertices[vertex_1];
    ListGraph::Node &y = lemon_vertices[vertex_2];

    // delete parallel edges of larger costs
    #ifdef DEBUG_CONTRACTION_WRAPPER
        cout << "deleting " << parallel_edges.size() << " edges:" << endl;
    #endif

    vector<ListGraph::Edge>::iterator parallel = parallel_edges.begin();
    while (parallel != parallel_edges.end())
    {
        #ifdef DEBUG_CONTRACTION_WRAPPER
            cout << "- {" << g.id(g.u(*parallel)) << "," << g.id(g.v(*parallel)) << "}, weight " << (*lemon_weight)[*parallel] << endl;
        #endif

        // index in the original graph
        dropped_indices.push_back( (*lemon_edges_inverted_index)[*parallel] );

        g.erase(*parallel); // NB! INVALIDATES NOT ONLY ITERATORS, BUT ALSO USING "PARALLEL" AS KEY IN MAPS

        ++parallel;
    }

    #ifdef DEBUG_CONTRACTION_WRAPPER
        cout << "done!" << endl;
    #endif

    /*
    cout << "after contracting {" << g.id(x) << "," << g.id(y) << "}: \t edges from " << g.id(x) << " are ";

    for (ListGraph::IncEdgeIt edge_from_x(g,x); edge_from_x != INVALID; ++edge_from_x)
    {
            cout << g.id(g.u(edge_from_x)) << " to ";
            cout << g.id(g.v(edge_from_x)) << ", ";
    }
    cout << " \t and edges from " << g.id(y) << " are ";
    for (ListGraph::IncEdgeIt edge_from_y(g,y); edge_from_y != INVALID; ++edge_from_y)
    {
            cout << g.id(g.u(edge_from_y)) << " to ";
            cout << g.id(g.v(edge_from_y)) << ", ";
    }
    cout << endl;
    */

    // TO DO: REDESIGN THESE OPERATIONS, AS THIS IS ERROR-PRONE, IN NEW APPLICATIONS! ONLY THE LEMON DATA-STRUCTURE IS CORRECT AFTERWARDS
    // vertex y is deleted next, so we fix the edge list so that edges previously incident to y now incide on x
    //cout << "contracting edge {" << vertex_1 << "," << vertex_2 << "}" << endl;
    for (ListGraph::IncEdgeIt edge_from_y(g,y); edge_from_y != INVALID; ++edge_from_y)
    {
        long edge_idx = (*lemon_edges_inverted_index)[edge_from_y];
        if (s[edge_idx] == vertex_2)
            s[edge_idx] = vertex_1;
        else
            t[edge_idx] = vertex_1;
        
        /*
        cout << "idx = " << edge_idx
             << ", s[" << edge_idx << "] = " << s[edge_idx]
             << ", t[" << edge_idx << "] = " << t[edge_idx]
             << endl << endl;
        */
    }

    g.contract(x, y, true);

    return dropped_indices;
}

bool Graph::mst()
{
    /***
     * Using LEMON to calculate a minimum spanning tree via their efficient 
     * implementation of Kruskal's algorithm (actually, the efficient union-find
     * is the game-changer). The MST cost and solution are stored in this object.
     */

    #ifdef DEBUG_MST
        cout << "Determining an MST in the graph with edges:" << endl;
        for (ListGraph::EdgeIt e_it(*lemon_graph); e_it != INVALID; ++e_it)
        {
            cout << "edge " << lemon_graph->id(e_it) << " is {"
                 << lemon_graph->id(lemon_graph->u(e_it)) << ","
                 << lemon_graph->id(lemon_graph->v(e_it)) << "}, ";
            cout << "inverted index = " << (*lemon_edges_inverted_index)[e_it] << ", ";
            cout << "weight = " << (*lemon_weight)[e_it] << endl;
        }
    #endif

    Timer mst_timer;
    mst_timer.start();
    vector<ListGraph::Edge> tree_edges;
    this->mst_weight = kruskal(*lemon_graph, *lemon_weight, back_inserter(tree_edges));

    // store 0-1 vector indicating which edges are in the MST
    this->mst_vector.clear();
    this->mst_vector = vector<bool>(num_edges, false);

    for (vector<ListGraph::Edge>::iterator it = tree_edges.begin(); it != tree_edges.end(); ++it)
    {
        long edge_index = (*lemon_edges_inverted_index)[*it];
        this->mst_vector[edge_index] = true;
    }

    mst_timer.halt();
    this->mst_runtime = mst_timer.realTime();

    #ifdef DEBUG_MST
        cout << "MST weight = " << this->mst_weight << endl;

        cout << "MST edges (in LEMON): " << endl;
        for (unsigned i = 0; i != tree_edges.size(); ++i)
            cout << "{" << lemon_graph->id(lemon_graph->u(tree_edges[i])) << ", "
                 << lemon_graph->id(lemon_graph->v(tree_edges[i])) << "}"
                 << ", Graph edge " << (*lemon_edges_inverted_index)[ tree_edges[i] ] << endl;

        cout << "MST vector: " << endl;
        for (unsigned i = 0; i != mst_vector.size(); ++i)
            cout << mst_vector[i];
        cout << endl;

        cout << "MST runtime: " << mst_timer.realTime() << endl;
    #endif

    if (tree_edges.size() != (unsigned) lemon::countNodes(*lemon_graph) - 1)
    {
        // KRUSKAL RETURNED A FOREST, NOT A TREE
        return false;
    }
    else
        return true;
}

bool Graph::maxst()
{
    /***
     * Using LEMON to calculate a maximum weight spanning tree (find MST using
     * negative/opposite edge weights). The corresponding weight is stored in
     * this object.
     */

    #ifdef DEBUG_MAXWEIGHTST
        cout << "Determining an MaxWeightST in the graph with edges:" << endl;
        for (ListGraph::EdgeIt e_it(*lemon_graph); e_it != INVALID; ++e_it)
        {
            cout << "edge " << lemon_graph->id(e_it) << " is {"
                 << lemon_graph->id(lemon_graph->u(e_it)) << ","
                 << lemon_graph->id(lemon_graph->v(e_it)) << "}, ";
            cout << "inverted index = " << (*lemon_edges_inverted_index)[e_it] << ", ";
            cout << "weight = " << (*opposite_weights)[e_it] << endl;
        }
    #endif

    vector<ListGraph::Edge> tree_edges;
    this->maxst_weight 
        = (-1)*kruskal(*lemon_graph, *opposite_weights, back_inserter(tree_edges));

    #ifdef DEBUG_MAXWEIGHTST
        cout << "MaxWeightST weight = " << this->maxst_weight << endl;

        cout << "MaxWeightST edges (in LEMON): " << endl;
        for (unsigned i = 0; i != tree_edges.size(); ++i)
            cout << "{" << lemon_graph->id(lemon_graph->u(tree_edges[i])) << ", "
                 << lemon_graph->id(lemon_graph->v(tree_edges[i])) << "}"
                 << ", Graph edge " << (*lemon_edges_inverted_index)[ tree_edges[i] ] << endl;
    #endif

    if (tree_edges.size() != (unsigned) lemon::countNodes(*lemon_graph) - 1)
    {
        // KRUSKAL RETURNED A FOREST, NOT A TREE
        return false;
    }
    else
        return true;
}

pair<bool,double> Graph::mst_probing_var(long probe_idx, bool probe_value)
{
    /***
     * Returns the cost of an MST removing an edge (if probe_value = 0) or
     * fixing it in the solution (if probe_value = 1). Since the former case
     * might result in a disconnected graph, the first element in the returned
     * pair indicates whether the bound in the second element corresponds to an
     * actual spanning tree (else, it corresponds to a minimum spanning forest)
     */

    #ifdef DEBUG_MST_PROBING
        cout << "probing MST with edge " << probe_idx;
        if (!probe_value) cout << " forbidden" << endl;
        else cout << " forced into the tree" << endl;

        cout << "original graph:" << endl;
        for (ListGraph::EdgeIt e_it(*lemon_graph); e_it != INVALID; ++e_it)
        {
            cout << "edge " << lemon_graph->id(e_it) << " is {"
                 << lemon_graph->id(lemon_graph->u(e_it)) << ","
                 << lemon_graph->id(lemon_graph->v(e_it)) << "}, ";
            cout << "inverted index = " << (*lemon_edges_inverted_index)[e_it] << ", ";
            cout << "mst_weight = " << (*lemon_weight)[e_it] << endl;
        }
    #endif

    vector<ListGraph::Edge> mst_edges;
    double mst_weight;
    Timer probing_timer;
    probing_timer.start();

    long vertex_1 = this->s[probe_idx];
    long vertex_2 = this->t[probe_idx];

    // Case 1: probing MST forcing the edge
    if (probe_value == true)
    {
        /* 1.1 DUPLICATE THE LEMON ADJACENCY LIST
         * Only this structure is edited next! Using the LEMON class template
         * GraphCopy to this end.
         */

        Timer copy_timer;
        copy_timer.start();

        ListGraph graph_copy;

        GraphCopy<ListGraph, ListGraph> copier(*this->lemon_graph, graph_copy);

        // save vertex correspondence: must be from old to new
        ListGraph::NodeMap<ListGraph::Node> vertex_ref_old2new(*this->lemon_graph);
        copier.nodeRef(vertex_ref_old2new);

        // save inverted edge correspondences: must be from new to old
        ListGraph::EdgeMap<ListGraph::Edge> edge_ref_new2old(graph_copy);
        copier.edgeCrossRef(edge_ref_new2old);

        // copy edge map
        ListGraph::EdgeMap<double> weights_copy(graph_copy);
        copier.edgeMap(*this->lemon_weight, weights_copy);

        copier.run();   // actually executes the copy actions determined above

        copy_timer.halt();

        #ifdef DEBUG_MST_PROBING
            cout << "copy runtime " << copy_timer.realTime() << endl;
        #endif

        // 1.2 CONTRACT EDGE
        vector<ListGraph::Edge> parallel_edges = 
        lemon_parallel_edges_if_contract( graph_copy,
                                          weights_copy,
                                          vertex_ref_old2new[lemon_vertices[vertex_1]],
                                          vertex_ref_old2new[lemon_vertices[vertex_2]] );

        // references to make notation less cumbersome
        ListGraph &g = graph_copy;
        ListGraph::Node &x = vertex_ref_old2new[lemon_vertices[vertex_1]];
        ListGraph::Node &y = vertex_ref_old2new[lemon_vertices[vertex_2]];

        // delete parallel edges of larger costs
        #ifdef DEBUG_CONTRACTION_WRAPPER
            cout << "deleting " << parallel_edges.size() << " edges:" << endl;
        #endif

        vector<ListGraph::Edge>::iterator parallel = parallel_edges.begin();
        while (parallel != parallel_edges.end())
        {
            #ifdef DEBUG_CONTRACTION_WRAPPER
                cout << "- {" << g.id(g.u(*parallel)) << "," << g.id(g.v(*parallel)) << "}, weight " << weights_copy[*parallel] << endl;
            #endif

            g.erase(*parallel); // NB! INVALIDATES NOT ONLY ITERATORS, BUT ALSO USING "PARALLEL" AS KEY IN MAPS

            ++parallel;
        }

        #ifdef DEBUG_CONTRACTION_WRAPPER
            cout << "done!" << endl;
        #endif

        g.contract(x, y, true);

        #ifdef DEBUG_MST_PROBING
            cout << "the edited graph copy consists of:" << endl;
            for (ListGraph::EdgeIt e_it(graph_copy); e_it != INVALID; ++e_it)
            {
                cout << "edge " << graph_copy.id(e_it) << " is {"
                     << graph_copy.id(graph_copy.u(e_it)) << ","
                     << graph_copy.id(graph_copy.v(e_it)) << "}, ";
                cout << "original inverted index = " << (*lemon_edges_inverted_index)[ edge_ref_new2old[e_it] ] << ", ";
                cout << "weight = " << weights_copy[e_it] << endl;
            }
        #endif

        // 1.3 RUN MST
        mst_weight = kruskal(graph_copy, weights_copy, back_inserter(mst_edges));
    }

    // Case 2: probing MST without the edge
    else
    {
        /* 2.1 DELETE EDGE IN THE ORIGINAL GRAPH (AVOID COPYING GRAPH, AND
         * TRAVERSING A NEIGHBOURHOOD TO RETRIEVE EDGE OBJECT IN THE COPY)
         */
        lemon_graph->erase(this->lemon_edges[probe_idx]);

        // 2.2 RUN MST
        mst_weight = kruskal(*this->lemon_graph, *this->lemon_weight, back_inserter(mst_edges));

        // 2.3 REINSERT EDGE
        ListGraph::Edge e = this->lemon_graph->addEdge(this->lemon_vertices[vertex_1], this->lemon_vertices[vertex_2]);
        this->lemon_edges[probe_idx] = e;
        (*this->lemon_weight)[e] = this->w[probe_idx];
        (*this->lemon_edges_inverted_index)[e] = probe_idx;
    }

    // 3. PREPARE RETURN VALUE
    bool result_is_tree = true;
    double result_weight = mst_weight;

    if (probe_value == true)
    {
        // add cost of contracted edge
        result_weight += this->w[probe_idx];
    }
    else
    {
        // check if solution after we deleted the edge is a spanning tree or forest
        if (mst_edges.size() != (unsigned) lemon::countNodes(*lemon_graph) - 1)
        {
            result_is_tree = false;
            #ifdef DEBUG_MST_PROBING
                cout << "WARNING: kruskal returns a spanning forest (" << mst_edges.size() << " edges!)" << endl;
            #endif
        }
    }

    #ifdef DEBUG_MST_PROBING
        cout << "original graph should be intact:" << endl;
        for (ListGraph::EdgeIt e_it(*lemon_graph); e_it != INVALID; ++e_it)
        {
            cout << "edge " << lemon_graph->id(e_it) << " is {"
                 << lemon_graph->id(lemon_graph->u(e_it)) << ","
                 << lemon_graph->id(lemon_graph->v(e_it)) << "}, ";
            cout << "inverted index = " << (*lemon_edges_inverted_index)[e_it] << ", ";
            cout << "mst_weight = " << (*lemon_weight)[e_it] << endl;
        }
    #endif

    // runtime saved on the object, for eventual use later on
    probing_timer.halt();
    this->probe_runtime = probing_timer.realTime();

    return make_pair(result_is_tree, result_weight);
}

ListGraph::Edge Graph::lemon_test_adj_getting_edge(ListGraph &g, ListGraph::Node &x, ListGraph::Node &y)
{
    /***
     * Auxiliary function to test adjacency in the LEMON data structure.
     * Returns edge, if found; otherwise, returns INVALID
     */

    for (ListGraph::IncEdgeIt e(g, x); e != INVALID; ++e)
    {
        if ( g.id(g.v(e)) == g.id(y) || g.id(g.u(e)) == g.id(y))
            return e;
    }
    return INVALID;
}

bool Graph::lemon_test_adj(ListGraph &g, ListGraph::Node &x, ListGraph::Node &y)
{
    /// auxiliary function to test adjacency in the LEMON data structure

    for (ListGraph::IncEdgeIt e(g, x); e != INVALID; ++e)
    {
        if ( g.id(g.v(e)) == g.id(y) || g.id(g.u(e)) == g.id(y))
            return true;
    }
    return false;
}

vector<ListGraph::Edge> 
Graph::lemon_parallel_edges_if_contract(ListGraph &g,
                                        ListGraph::EdgeMap<double> &w,
                                        ListGraph::Node &x,
                                        ListGraph::Node &y)
{
    /**
     * wrapper over LEMON, contracting edge xy and keeping only the
     * least-weight (given by w) of any parallel edges
     */

    vector<ListGraph::Edge> deleted_edges;
    
    #ifdef DEBUG_CONTRACTION_WRAPPER
        cout << "contracting edge {" << g.id(x) << "," << g.id(y) << "}" << endl;
    #endif
    
    // traverses both neighbourhoods searching for common neighbours
    for (ListGraph::IncEdgeIt edge_from_x(g,x); edge_from_x != INVALID; ++edge_from_x)
    {
        // other end of this edge is not y
        if ( g.id(g.u(edge_from_x)) != g.id(y) && g.id(g.v(edge_from_x)) != g.id(y) )
        {
            // let z be the other end of this edge {x,z}
            #ifdef DEBUG_CONTRACTION_WRAPPER
                cout << "edge {" << g.id(g.u(edge_from_x)) << ","
                     << g.id(g.v(edge_from_x)) << "}, of weight " << w[edge_from_x]
                     << "... so z=";
            #endif
            
            ListGraph::Node z = g.id(g.u(edge_from_x)) != g.id(x) ?
                                g.u(edge_from_x) : g.v(edge_from_x) ;
            
            #ifdef DEBUG_CONTRACTION_WRAPPER
                cout << g.id(z) << " - ";
            #endif

            // search for z in N(y)
            bool z_neighbour_of_y = false;
            ListGraph::IncEdgeIt edge_from_y(g,y);
            while (edge_from_y != INVALID && !z_neighbour_of_y)
            {
                // let w be the other end of this edge {y,w}
                ListGraph::Node w = g.id(g.u(edge_from_y)) != g.id(y) ?
                                    g.u(edge_from_y) : g.v(edge_from_y) ;
                if ( g.id(z) == g.id(w) )
                    z_neighbour_of_y = true;
                else
                    ++edge_from_y;
            }

            if (z_neighbour_of_y)
            {
                #ifdef DEBUG_CONTRACTION_WRAPPER
                    cout << "also a neighbour from y!" << endl;
                    cout << "edge {" << g.id(g.u(edge_from_y)) << "," 
                         << g.id(g.v(edge_from_y)) << "}, of weight "
                         << w[edge_from_y] << "... ";
                #endif
                
                if (w[edge_from_y] >= w[edge_from_x])
                {
                    #ifdef DEBUG_CONTRACTION_WRAPPER
                        cout << "deleted!" << endl;
                    #endif
                    
                    deleted_edges.push_back(edge_from_y);
                }
                else
                {
                    #ifdef DEBUG_CONTRACTION_WRAPPER
                        cout << "kept!" << endl;
                    #endif

                    deleted_edges.push_back(edge_from_x);
                }
            }
            else
            {
                #ifdef DEBUG_CONTRACTION_WRAPPER
                    cout << "not a neighour from y" << endl;
                #endif
            }
        }
    }

    return deleted_edges;
}
