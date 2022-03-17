#include "io.h"

IO::IO()
{
    // only for consistency (if IO is used without calling parse_input_file)
    this->graph = new Graph();

    this->num_conflicts = 0;
}

IO::~IO()
{
    delete graph;
    conflicts.clear();
    conflict_graph_adj_list.clear();
}

// parse instance file
bool IO::parse_input_file(string filename)
{
    // check if this is an instance in the original benchmark (gcc extension)
    size_t dot_pos = filename.find_last_of(".");
    bool gcc_format = filename.find("gcc", dot_pos+1) != string::npos;

    ifstream input_fh(filename);
    
    if (input_fh.is_open())
    {
        // in gcc files, the first lines might have comments, followed by an instance id
        if (gcc_format)
        {
            string line;
            
            // skip comment lines
            do
            {
                getline(input_fh, line);
            }
            while(line.find("#") != string::npos);   // there is a trail

            instance_id.assign(line);
        }
        else
            instance_id.assign(filename);

        // trimmed instance id: contents after last slash and before the last dot
        size_t last_slash_pos = filename.find_last_of("/\\");
        instance_id_trimmed = filename.substr(last_slash_pos+1, dot_pos-1 - last_slash_pos);

        // 3 lines for cardinalities of vertex, edge and conflict sets
        long num_vertices, num_edges;
        input_fh >> num_vertices;
        input_fh >> num_edges;
        input_fh >> this->num_conflicts;

        // initialize graph (own structures and lemon object)
        delete graph;
        this->graph = new Graph(num_vertices,num_edges);
        this->graph->init_index_matrix();
        this->graph->init_lemon();

        // conflict graph structures
        this->conflicts.reserve(this->num_conflicts);
        this->conflict_graph_adj_list.reserve(num_edges);
        this->conflict_graph_adj_list.insert(this->conflict_graph_adj_list.begin(), 
            num_edges, list<long>() );   // adds num_edges copies of an empty list
        
        // m lines for edges in the instance graph
        for (long line_idx=0; line_idx<num_edges; ++line_idx)
        {
            // read edge terminal vertices and weight: i j w  (incidentally, i<j)
            long i, j;
            double w;
            
            input_fh >> i;
            graph->s.push_back(i);

            input_fh >> j;
            graph->t.push_back(j);
            
            input_fh >> w;
            graph->w.push_back(w);

            graph->adj_list[i].push_back(j);
            graph->adj_list[j].push_back(i);
            
            // should never happen
            if (graph->index_matrix[i][j] >= 0 ||
                graph->index_matrix[j][i] >= 0 )
            {
                cerr << "ERROR: repeated edge in input file line " << line_idx << endl << endl;
                return false;
            }
            
            // store index of current edge
            graph->index_matrix[i][j] = line_idx;
            graph->index_matrix[j][i] = line_idx;

            // lemon edge
            ListGraph::Edge e = graph->lemon_graph->addEdge(graph->lemon_vertices[i], graph->lemon_vertices[j]);
            graph->lemon_edges.push_back(e);
            (*graph->lemon_weight)[e] = w;
            (*graph->opposite_weights)[e] = (-1)*w;
            (*graph->lemon_edges_inverted_index)[e] = line_idx;
        }

        /* original instances (gcc format) include a single line per conflict,
         * but carrabs instances (cms format) include two lines for each (both
         * e1 e2 and e2 e1)
         */

        long conflict_lines_to_read = gcc_format ? num_conflicts : 2*num_conflicts;

        // p lines for conflicting edge pairs in the instance
        for (long line_idx=0; line_idx<conflict_lines_to_read; ++line_idx)
        {
            /* read terminal nodes from each edge in the pair:
             * a b c d (incidentally, a<b and c<d)
             */
            long a, b, c, d;
            input_fh >> a;
            input_fh >> b;
            input_fh >> c;
            input_fh >> d;
            
            long e1_index = graph->index_matrix[a][b];
            long e2_index = graph->index_matrix[c][d];

            // should never happen
            if (e1_index<0 || e2_index<0)
            {
                cerr << "ERROR: edge not found in conflict line " << line_idx << endl << endl;
                return false;
            }
            
            bool repeated = false;
            list<long>::iterator neighbour_it = this->conflict_graph_adj_list[e1_index].begin();
            while ( !repeated && neighbour_it != this->conflict_graph_adj_list[e1_index].end() )
            {
                if (*neighbour_it == e2_index)
                    repeated = true;

                ++neighbour_it;
            }

            // store conflict
            if (!repeated)
            {
                pair<long,long> edges = (a<c) ? make_pair(e1_index, e2_index) : make_pair(e2_index, e1_index);

                this->conflicts.push_back(edges);

                // update adjacency lists
                this->conflict_graph_adj_list[e1_index].push_back(e2_index);
                this->conflict_graph_adj_list[e2_index].push_back(e1_index);
            }
        }

        if ( (unsigned) num_conflicts != conflicts.size())
            cout << "num_conflicts = " << num_conflicts << ", conflict.size() = " << conflicts.size() << endl;
        
        input_fh.close();
    }
    else
    {
        cerr << "ERROR: could not open file (might not exist)." << endl;
        return false;
    }


    /*
    #ifdef DEBUG
        cout << "~~~ LEMON says:  " << endl;
        cout << "n = " << countNodes(*graph->lemon_graph) << endl;
        cout << "m = " << countEdges(*graph->lemon_graph) << endl;

        // iterating over edges
        for (ListGraph::EdgeIt e_it(*graph->lemon_graph); e_it != INVALID; ++e_it)
        {
            cout << "edge " << graph->lemon_graph->id(e_it) << " is {"
                 << graph->lemon_graph->id(graph->lemon_graph->u(e_it))
                 << "," << graph->lemon_graph->id(graph->lemon_graph->v(e_it)) << "}, ";
            cout << "inverted index = " << (*graph->lemon_edges_inverted_index)[e_it] << ", ";
            cout << "weight = " << (*graph->lemon_weight)[e_it] << endl;
        }

        // iterating over vertices and their neighbourhood
        for (ListGraph::NodeIt vertex(*graph->lemon_graph); vertex != INVALID; ++vertex)
        {
            int cnt = 0;
            for (ListGraph::IncEdgeIt e_it(*graph->lemon_graph,vertex); e_it != INVALID; ++e_it)
                cnt++;

            cout << "deg(" << graph->lemon_graph->id(vertex) << ") = " << cnt << endl;
        }
        cout << "~~~ ." << endl;
    #endif
    */

    return true;
}

bool IO::test_stability(vector<bool> &point)
{
    /// tests if an incidence vector satisfies all conflict constraints
    for ( vector< pair<long,long> >::iterator it = this->conflicts.begin();
          it != this->conflicts.end(); ++it )
    {
        long e1 = (*it).first;
        long e2 = (*it).second;
        if (point[e1] && point[e2])
            return false;
    }

    return true;
}

bool IO::test_acyclic(vector<bool> &point)
{
    /// tests if an incidence vector induces an acyclic subgraph

    bool acyclic = true;

    vector<bool> check = vector<bool>(graph->num_vertices, false);
    long checked_count = 0;

    // the search starts only once, if the subgraph is connected/acyclic
    long source = 0;
    while (source < graph->num_vertices && acyclic && checked_count < graph->num_vertices)
    {
        #ifdef DEBUG_DFS
            cout << "source = " << source << endl;
        #endif

        if (!check[source])
        {
            #ifdef DEBUG_DFS
                cout << "dfs starting!" << endl << endl;
            #endif
            dfs_checking_acyclic(source, -1, check, checked_count, point, acyclic);
        }

        ++source;
    }

    #ifdef DEBUG_DFS
        cout << "acyclic = " << acyclic << endl;
        cout << "check vector: { ";
        for (long i=0; i<graph->num_vertices; ++i) cout << check[i] << " ";
        cout << "}" << endl;
        cout << "checked_count = " << checked_count << endl;
    #endif

    return acyclic;
}

bool IO::test_acyclic_kstab(vector<bool> &point)
{
    /***
     * Tests if an incidence vector of num_vertices - 1 edges induces an
     * acyclic/connected subgraph.
     */

    bool acyclic = true;

    vector<bool> check = vector<bool>(graph->num_vertices, false);
    long checked_count = 0;

    // the search starts only once, if the subgraph is connected/acyclic
    dfs_checking_acyclic(0, -1, check, checked_count, point, acyclic);

    #ifdef DEBUG_DFS
        cout << "acyclic = " << (checked_count == graph->num_vertices) << endl;
        cout << "check vector: { ";
        for (long i=0; i<graph->num_vertices; ++i) cout << check[i] << " ";
        cout << "}" << endl;
        cout << "checked_count = " << checked_count << endl;
    #endif

    return (checked_count == graph->num_vertices);
}

void IO::dfs_checking_acyclic(long u,
                              long parent,
                              vector<bool> &check,
                              long &checked_count,
                              vector<bool> &active_edges,
                              bool &acyclic)
{
    /// dfs to check if the subgraph induced by active edges is acyclic

    check[u] = true;
    ++checked_count;

    list<long>::iterator it = graph->adj_list[u].begin();
    while (it != graph->adj_list[u].end() && acyclic)
    {
        long v = (*it);
        
        // considering only active edges
        long idx = graph->index_matrix[u][v];
        if (active_edges[idx])
        {
            if (!check[v])
                dfs_checking_acyclic(v, u, check, checked_count, active_edges, acyclic);
            else if (v != parent)
            {
                // vertex (not the antecessor) already visited: back edge!
                acyclic = false;
            }
        }

        ++it;
    }
}
