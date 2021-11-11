#include "io.h"

IO::IO()
{
    // only for consistency (if IO is used without calling parse_gcclib)
    this->graph = new Graph();

    this->num_conflicts = 0;
}

IO::~IO()
{
    delete graph;
    conflicts.clear();
    conflict_graph_adj_list.clear();
}

// parse GCCLib instance file
bool IO::parse_gcclib(const char *filename)
{
    ifstream input_fh(filename);
    
    if (input_fh.is_open())
    {
        string line;
        
        // skip comment lines
        do
        {
            getline(input_fh, line);
        }
        while(line.find("#") != string::npos);   // there is a trail
        
        // 1st line: instance id
        instance_id.assign(line);
        
        // 3 lines for cardinalities of vertex, edge and conflict sets
        long num_vertices, num_edges;
        input_fh >> num_vertices;
        input_fh >> num_edges;
        input_fh >> this->num_conflicts;

        // initialize graph adjacency list, edge list and lemon object
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
            // read edge terminal nodes and weight: i j w  (incidentally, i<j)
            long i, j, w;
            
            input_fh >> i;
            graph->s.push_back(i);

            input_fh >> j;
            graph->t.push_back(j);
            
            input_fh >> w;
            graph->w.push_back(w);
            
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
            (*graph->lemon_edges_inverted_index)[e] = line_idx;
        }
        
        // p lines for conflicting edge pairs in the instance
        for (long line_idx=0; line_idx<num_conflicts; ++line_idx)
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
            
            // store conflict
            pair<long,long> edges = (a<c) ? make_pair(e1_index, e2_index) : make_pair(e2_index, e1_index);
            
            // should never happen
            for (unsigned k=0; k<conflicts.size(); ++k)
            {
                pair<long, long> cur = conflicts[k];
                if (cur.first == edges.first || cur.first == edges.second)
                {
                    if (cur.second == edges.first || cur.second == edges.second)
                    {
                        cerr << "ERROR: repeated conflict in input file line " << line_idx << endl << endl;
                        return false;
                    }
                }
            }
            
            this->conflicts.push_back(edges);

            // update adjacency lists
            this->conflict_graph_adj_list[e1_index].push_back(e2_index);
            this->conflict_graph_adj_list[e2_index].push_back(e1_index);
        }
        
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
