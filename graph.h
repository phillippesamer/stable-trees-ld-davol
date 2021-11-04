#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>
#include <iostream>

#include <lemon/list_graph.h>
#include <lemon/kruskal.h>
#include <lemon/time_measure.h>
#include <lemon/core.h>

using namespace std;
using namespace lemon;

/***
 * \file graph.h
 * 
 * Module containing different data structures representing the same graph:
 * - an edge list, with either terminals of each edge stored in vector s, t
 * - an adjacency matrix storing the corresponding index (or -1 for non-edges)
 * - an adjacency list from the LEMON (Library for Efficient Modeling and 
 * Optimization in Networks), so as to use the highly efficient implementations
 * of algorithms they offer
 * 
 * IO and Model classes are declared friends to avoid cumbersome get/set methods.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 02.11.2021
 */
class Graph
{
public:
    Graph();
    Graph(long, long);
    ~Graph();

    void init_index_matrix();
    void free_index_matrix();

    void init_lemon();

private:
    friend class IO;
    friend class Model;

    long num_vertices;
    long num_edges;

    // edge list
    vector<long> s;        // terminal node 1
    vector<long> t;        // terminal node 2
    vector<long> w;        // edge weight

    // adjacency matrix storing edge indexes
    bool using_matrix;
    long **index_matrix;

    // LEMON adjacency list: http://lemon.cs.elte.hu/pub/doc/1.3/a00237.html
    bool using_lemon;
    ListGraph *lemon_graph;
    vector<ListGraph::Node> lemon_vertices;
    vector<ListGraph::Edge> lemon_edges;
    ListGraph::EdgeMap<long> *lemon_weight;
};

#endif
