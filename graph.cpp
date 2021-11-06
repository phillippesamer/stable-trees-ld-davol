#include "graph.h"

Graph::Graph()
{
    using_matrix = false;
    using_lemon = false;
    num_vertices = 0;
    num_edges = 0;
}

Graph::Graph(long n, long m)
{
    using_matrix = false;
    using_lemon = false;
    num_vertices = n;
    num_edges = m;

    this->s.reserve(m);
    this->t.reserve(m);
    this->w.reserve(m);
    this->mst_edges.reserve(n); // n-1, actually...
}

Graph::~Graph()
{
    s.clear();
    t.clear();
    w.clear();
    mst_edges.clear();

	if (using_matrix)
        free_index_matrix();

    if (using_lemon)
    {
        lemon_vertices.clear();
        lemon_edges.clear();
        delete lemon_weight;
        delete lemon_inverted_edge_index;
        delete lemon_graph;
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
    lemon_weight = new ListGraph::EdgeMap<long>(*lemon_graph);
    lemon_inverted_edge_index = new ListGraph::EdgeMap<long>(*lemon_graph);
}

void Graph::update_single_weight(long idx, long new_weight)
{
    // weight in the edge list
    this->w[idx] = new_weight;

    // weight in lemon's adjacency list
    ListGraph::Edge e = this->lemon_edges[idx];
    (*lemon_weight)[e] = new_weight;
}

void Graph::update_all_weights(vector<long> new_weights)
{
    // weight in the edge list
    this->w.clear();
    this->w = vector<long>(new_weights);

    // weight in lemon's adjacency list
    for (long i=0; i<num_edges; ++i)
    {
        ListGraph::Edge e = this->lemon_edges[i];
        (*lemon_weight)[e] = new_weights[i];
    }
}

bool Graph::mst()
{
    /***
     * Using LEMON to calculate a minimum spanning tree via their efficient 
     * implementation of Kruskal's algorithm (actually, the efficient union-find
     * is the game-changer).
     * The MST cost and solution are store in this object.
     */

    #ifdef DEBUG
        cout << "Determining an MST in the graph with edges:" << endl;
        for (ListGraph::EdgeIt e_it(*lemon_graph); e_it != INVALID; ++e_it)
        {
            cout << "edge " << lemon_graph->id(e_it) << " is {" << lemon_graph->id(lemon_graph->u(e_it)) << "," << lemon_graph->id(lemon_graph->v(e_it)) << "}, ";
            cout << "inverted index = " << (*lemon_inverted_edge_index)[e_it] << ", ";
            cout << "weight = " << (*lemon_weight)[e_it] << endl;
        }
    #endif

    Timer mst_timer;
    mst_timer.start();
    vector<ListGraph::Edge> tree_edges;
    this->mst_weight = kruskal(*lemon_graph, *lemon_weight, back_inserter(tree_edges));
    mst_timer.halt();

    #ifdef DEBUG
        cout << "MST weight = " << this->mst_weight << endl;
        cout << "MST edges: " << endl;
        for (unsigned i = 0; i != tree_edges.size(); ++i)
            cout << "{" << lemon_graph->id(lemon_graph->u(tree_edges[i])) << ", " << lemon_graph->id(lemon_graph->v(tree_edges[i])) << "}" << endl;
        cout << "Elapsed time: " << mst_timer.realTime() << endl;
    #endif

    // store MST edges from lemon_inverted_edge_index
    this->mst_edges.clear();
    for (vector<ListGraph::Edge>::iterator it = tree_edges.begin(); it != tree_edges.end(); ++it)
            this->mst_edges.push_back( (*lemon_inverted_edge_index)[*it] );

    return true;
}

bool Graph::lemon_test_adj(ListGraph &g, ListGraph::Node &x, ListGraph::Node &y)
{
    /// auxiliary function to test adjacency in the LEMON data structure

    for (ListGraph::IncEdgeIt e(g, x); e != INVALID; ++e)
    {
        //cout << "testing if " << g.id(y) << " is equal to " << g.id(g.v(e)) << " or " << g.id(g.u(e)) << endl;
        if ( g.id(g.v(e)) == g.id(y) || g.id(g.u(e)) == g.id(y))
            return true;
    }
    return false;
}

long lemon_contract_dropping_parallel_edges(ListGraph &g, ListGraph::EdgeMap<long> &w, ListGraph::Node &x, ListGraph::Node &y)
{
    /**
     * wrapper over LEMON, contracting edge xy and keeping only the
     * least-weight (given by w) of any parallel edges
     */

    long deleted_edges = 0;
    vector<ListGraph::Edge> to_be_deleted;
    
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
            cout << "edge {" << g.id(g.u(edge_from_x)) << "," << g.id(g.v(edge_from_x)) << "}, of weight " << w[edge_from_x] << "... so z=";
            #endif
            
            ListGraph::Node z = g.id(g.u(edge_from_x)) != g.id(x) ? g.u(edge_from_x) : g.v(edge_from_x) ;
            
            #ifdef DEBUG_CONTRACTION_WRAPPER
            cout << g.id(z) << " - ";
            #endif

            // search for z in N(y)
            bool z_neighbour_of_y = false;
            ListGraph::IncEdgeIt edge_from_y(g,y);
            while (edge_from_y != INVALID && !z_neighbour_of_y)
            {
                // let w be the other end of this edge {y,w}
                ListGraph::Node w = g.id(g.u(edge_from_y)) != g.id(y) ? g.u(edge_from_y) : g.v(edge_from_y) ;
                if ( g.id(z) == g.id(w) )
                    z_neighbour_of_y = true;
                else
                    ++edge_from_y;
            }

            if (z_neighbour_of_y)
            {
                #ifdef DEBUG_CONTRACTION_WRAPPER
                cout << "also a neighbour from y!" << endl;
                cout << "edge {" << g.id(g.u(edge_from_y)) << "," << g.id(g.v(edge_from_y)) << "}, of weight " << w[edge_from_y] << "... ";
                #endif
                
                if (w[edge_from_y] >= w[edge_from_x])
                {
                    #ifdef DEBUG_CONTRACTION_WRAPPER
                    cout << "deleted!" << endl;
                    #endif
                    
                    to_be_deleted.push_back(edge_from_y);
                }
                else
                {
                    #ifdef DEBUG_CONTRACTION_WRAPPER
                    cout << "kept!" << endl;
                    #endif

                    to_be_deleted.push_back(edge_from_x);
                }
                ++deleted_edges;
            }
            else
            {
                #ifdef DEBUG_CONTRACTION_WRAPPER
                cout << "not a neighour from y" << endl;
                #endif
            }
        }
    }

    // delete parallel edges of larger costs
    #ifdef DEBUG_CONTRACTION_WRAPPER
    cout << "deleting " << deleted_edges << "=" << to_be_deleted.size() << " edges:" << endl;
    #endif

    for (long i=0; i<deleted_edges; ++i)
    {
        #ifdef DEBUG_CONTRACTION_WRAPPER
        cout << "- {" << g.id(g.u(to_be_deleted[i])) << "," << g.id(g.v(to_be_deleted[i])) << "}, weight " << w[to_be_deleted[i]] << endl;
        #endif

        g.erase(to_be_deleted[i]);
    }

    #ifdef DEBUG_CONTRACTION_WRAPPER
    cout << "done!" << endl;
    #endif

    g.contract(x, y, true);

    return deleted_edges;
}
