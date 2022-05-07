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
    
    bool parse_input_file(string, bool);

    // record of effects due to edges fixed in preprocessing
    bool problem_solved;
    bool problem_infeasible;
    double objective_offset;
    vector<long> fixed_vars;
    vector<long> preprocessed_to_original_edge_idx;
    double preprocessing_runtime;

    bool run_mst() { return graph->mst(); }
    double get_mst_weight() { return graph->mst_weight; }
    double get_mst_runtime() { return graph->mst_runtime; }

    bool run_maxst() { return graph->maxst(); }
    double get_maxst_weight() { return graph->maxst_weight; }

    bool test_stability(vector<bool> &);
    bool test_acyclic(vector<bool> &);
    bool test_acyclic_kstab(vector<bool> &);

    // instance data
    long num_conflicts;
    string instance_id;
    string instance_id_trimmed;

private:
    friend class KStabModel;
    friend class StableSpanningTreeModel;
    friend class LDDA;
    friend class LDDAVolume;
    friend class KStabCutGenerator;
    friend class SSTCutGenerator;

    Graph *graph;  // different representations of the original graph

    vector< pair<long,long> > conflicts;  // conflicting edges (indices)

    vector< list<long> > conflict_graph_adj_list;

    void dfs_checking_acyclic(long, long, vector<bool>&, long&, vector<bool>&, bool&);

    // preprocessing algorithm methods
    bool preprocess();
    long preprocess_bridges();
    long preprocess_chain();
    long preprocess_pairwise();
    void remove_edge(long);
    void add_conflict(long, long);
    void dfs_for_bridges(long, long&, long, vector<long>&, bool*, long*, long*, vector< vector<long> >&);
    void dfs_for_chain(long, long&, long, vector<long>&, long*, long*, vector< vector<long> >&, bool*, bool*);

    void start_timer();
    double total_time();
    struct timeval* preprocessing_clock_start;
    struct timeval* preprocessing_clock_stop;
};

#endif
