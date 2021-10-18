#ifndef _IO_H_
#define _IO_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

using namespace std;

// simple edge list representation of the graph augmented with the conflicts
class IO
{
public:
    ~IO();
    
    bool parse_gcclib(const char *);
    
    // instance data
    int num_vertices;
    int num_edges;
    int num_conflicts;
    string instance_id;
    
    vector<int> s;        // terminal node 1
    vector<int> t;        // terminal node 2
    vector<int> w;        // edge weight
    
    int **index_matrix;   // adjacency matrix storing edge indexes
    
    vector< pair<int,int> > conflicts;   // conflicting edges (indexes only)
    
private:
    void init_index_matrix();
    void free_index_matrix();
};

#endif
