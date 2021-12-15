#ifndef _MSTCCMODEL_H_
#define _MSTCCMODEL_H_

#include <iostream>
#include <sstream>
#include <stack>
#include <sys/time.h>

#include <set>
#include <iomanip>

#include "gurobi_c++.h"

// using the preflow-push algorithm in LEMON (COIN-OR project, see: www.lemon.cs.elte.hu/)
#include <lemon/concepts/digraph.h>
#include <lemon/smart_graph.h>
#include <lemon/preflow.h>

using namespace lemon;
using namespace lemon::concepts;

#include "io.h"
#include "kstab_model.h"

#define VIOLATION_TOL 0.1

/***
 * \file mstcc_model.h
 * 
 * Module for the natural integer programming formulation to find stable
 * spanning trees of minimum weight, using the Gurobi solver API.
 * 
 * Extends the fixed cardinality stable set model class. Used in LDDA project
 * only to determine the LP relaxation bound.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 04.12.2021
 */

typedef struct
{
    map<string,int> pool;               // stores cuts not to add repeated ones

    vector<int> sec_diff_cuts;          // # of different cuts in each separation
    vector<double> sec_infeas_std_dev;  // std.dev. of the amount violated    
} statistics;


class StableSpanningTreeModel: public KStabModel
{
public:
    StableSpanningTreeModel(IO*);
    ~StableSpanningTreeModel();
    
    bool solve_lp_relax(bool); // full model (with SEC in the ORIGINAL GRAPH)
    double lp_bound;
    double lp_runtime;
    long lp_passes;

private:
    friend class violated_sec;   // encapsulates a class of violated cuts

    statistics *stats;

    bool separate_SEC(vector<GRBLinExpr> &, vector<long> &);
    bool separate_SEC_integer(vector<GRBLinExpr> &, vector<long> &);
    bool separate_SEC_fallback(vector<GRBLinExpr> &, vector<long> &);
    void dfs_checking_acyclic(long, long, vector<bool> &, long &, stack<long> &, vector< list<long> > &, bool &);
    void dfs_to_count(long, bool *, bool **, vector<long> &, vector< list<long> > &, long*, long*);
};

#endif
