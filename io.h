#ifndef _IO_H_
#define _IO_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <algorithm>

#include "graph.h"

using namespace std;

/***
 * \file io.h
 * 
 * Module for input and output functionality, including a Graph object for
 * main data structures.
 * 
 * Some classes are declared friends to avoid cumbersome get/set calls.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 01.11.2021
 */
class IO
{
public:
    IO();
    virtual ~IO();
    
    bool parse_gcclib(const char *);

    long num_vertices() { return graph->num_vertices; }
    long num_edges() { return graph->num_edges; }

    bool test_stability(vector<bool> &);
    bool test_acyclic(vector<bool> &);
    bool test_acyclic_kstab(vector<bool> &);

private:
    friend class KStabModel;
    friend class StableSpanningTreeModel;
    friend class LDDA;

    // instance data
    long num_conflicts;
    string instance_id;

    Graph *graph;  // different representations of the original graph

    vector< pair<long,long> > conflicts;  // conflicting edges (indices)

    vector< list<long> > conflict_graph_adj_list;

    void dfs_checking_acyclic(long, long, vector<bool>&, long&, vector<bool>&, bool&);
};

#endif
