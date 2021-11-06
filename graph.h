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
 * IO, Model and LDDA classes are declared friends to avoid cumbersome get/set calls.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 01.11.2021
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

    void update_single_weight(long,long);
    void update_all_weights(vector<long>);

    bool mst();
    long mst_weight;
    vector<long> mst_edges;

private:
    friend class IO;
    friend class Model;
    friend class LDDA;

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
    ListGraph::EdgeMap<long> *lemon_inverted_edge_index;
    ListGraph::EdgeMap<long> *lemon_weight;

    bool lemon_test_adj(ListGraph &, ListGraph::Node &, ListGraph::Node &);
    long lemon_contract_dropping_parallel_edges(ListGraph &, ListGraph::EdgeMap<long> &, ListGraph::Node &, ListGraph::Node &);
};

#endif
