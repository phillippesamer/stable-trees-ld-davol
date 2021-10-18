#include "io.h"
#include <algorithm>

IO::~IO()
{
    s.clear();
    t.clear();
    w.clear();
    conflicts.clear();
    
    free_index_matrix();
}

// adjacency matrix storing edge indexes
void IO::init_index_matrix()
{
    index_matrix = new int*[num_vertices];
    for (int i=0; i<num_vertices; ++i)
    {
        index_matrix[i] = new int[num_vertices];
        for (int j=0; j<num_vertices; ++j)
            index_matrix[i][j] = -1;
    }
}

void IO::free_index_matrix()
{
    for (int i = 0; i<num_vertices; ++i)
        delete[] index_matrix[i];

    delete[] index_matrix;
}

// parse GCCLib instance file into IO object
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
        input_fh >> num_vertices;
        input_fh >> num_edges;
        input_fh >> num_conflicts;
        
        // initialize index matrix in edge list
        init_index_matrix();
        
        // m lines for edges in the instance graph
        for (int line=0; line<num_edges; ++line)
        {
            // read edge terminal nodes and weight: i j w , such that i<j
            int i, j, w;
            
            input_fh >> i;
            this->s.push_back(i);

            input_fh >> j;
            this->t.push_back(j);
            
            input_fh >> w;
            this->w.push_back(w);
            
            // should never happen
            if (index_matrix[i][j] >= 0 ||
                index_matrix[j][i] >= 0 )
            {
                cerr << "ERROR: repeated edge in input file line " << line << endl << endl;
                return false;
            }
            
            // store index of current edge
            index_matrix[i][j] = line;
            index_matrix[j][i] = line;
        }
        
        // p lines for conflicting edge pairs in the instance
        for (int line=0; line<num_conflicts; ++line)
        {
            /* read terminal nodes from each edge in the pair: a b c d ,
             * such that a<b and c<d
             */
            int a, b, c, d;
            input_fh >> a;
            input_fh >> b;
            input_fh >> c;
            input_fh >> d;
            
            int e1_index = index_matrix[a][b];
            int e2_index = index_matrix[c][d];
            
            // should never happen
            if (e1_index<0 || e2_index<0)
            {
                cerr << "ERROR: edge not found in conflict line " << line << endl << endl;
                return false;
            }
            
            // store conflict
            pair<int,int> edges = (a<c) ? make_pair(e1_index, e2_index) : make_pair(e2_index, e1_index);
            
            // should never happen
            for (unsigned k=0; k<conflicts.size(); ++k)
            {
                pair<int, int> cur = conflicts[k];
                if (cur.first == edges.first || cur.first == edges.second)
                {
                    if (cur.second == edges.first || cur.second == edges.second)
                    {
                        cerr << "ERROR: repeated conflict in input file line " << line << endl << endl;
                        return false;
                    }
                }
            }
            
            conflicts.push_back(edges);
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
