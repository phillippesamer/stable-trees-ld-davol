#ifndef _MSTCCMODEL_H_
#define _MSTCCMODEL_H_

#include <iostream>
#include <sstream>
#include <iomanip>
#include <stack>
#include <set>
#include <sys/time.h>
#include <utility>

#include "gurobi_c++.h"

// using the preflow-push algorithm in LEMON (COIN-OR project, see: www.lemon.cs.elte.hu/)
#include <lemon/concepts/digraph.h>
#include <lemon/smart_graph.h>
#include <lemon/preflow.h>

using namespace lemon;
using namespace lemon::concepts;

#include "io.h"
#include "kstab_model.h"

#define SEC_VIOLATION_TOL 0.0001
#define SEC_SEPARATION_PRECISION 14   // <= 14 without changing everything to long double

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
    virtual ~StableSpanningTreeModel();
    
    bool solve_lp_relax(bool); // full model (with SEC in the ORIGINAL GRAPH)
    double lp_bound;
    double lp_runtime;
    long lp_passes;

private:
    friend class violated_sec;   // encapsulates a class of violated cuts

    statistics *stats;

    bool separate_SEC_integer(vector<GRBLinExpr> &, vector<long> &);
    bool separate_SEC_folklore(vector<GRBLinExpr> &, vector<long> &);
    bool separate_SEC_classical(vector<GRBLinExpr> &, vector<long> &);
    void dfs_checking_acyclic(long, long, vector<bool> &, long &, stack<long> &, vector< list<long> > &, bool &);
};

#endif
