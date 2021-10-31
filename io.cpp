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

        this->conflicts.reserve(this->num_conflicts);
        
        // initialize graph adjacency list, edge list and lemon object
        delete graph;
        this->graph = new Graph(num_vertices,num_edges);
        this->graph->init_index_matrix();
        
        // m lines for edges in the instance graph
        for (long line=0; line<num_edges; ++line)
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
                cerr << "ERROR: repeated edge in input file line " << line << endl << endl;
                return false;
            }
            
            // store index of current edge
            graph->index_matrix[i][j] = line;
            graph->index_matrix[j][i] = line;
        }
        
        // p lines for conflicting edge pairs in the instance
        for (long line=0; line<num_conflicts; ++line)
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
                cerr << "ERROR: edge not found in conflict line " << line << endl << endl;
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
                        cerr << "ERROR: repeated conflict in input file line " << line << endl << endl;
                        return false;
                    }
                }
            }
            
            this->conflicts.push_back(edges);
        }
        
        input_fh.close();
    }
    else
    {
        cerr << "ERROR: could not open file (might not exist)." << endl;
        return false;
    }
    
    return true;
}
