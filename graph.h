#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>

using namespace std;

class Graph
{
public:
    Graph();
    Graph(long, long);
    ~Graph();

    void init_index_matrix();
    void free_index_matrix();

private:
    friend class IO;
    friend class Model;

    long num_vertices;
    long num_edges;

    // edge list
    vector<long> s;        // terminal node 1
    vector<long> t;        // terminal node 2
    vector<long> w;        // edge weight
    
    long **index_matrix;   // adjacency matrix storing edge indexes
    bool using_matrix;
};

#endif
