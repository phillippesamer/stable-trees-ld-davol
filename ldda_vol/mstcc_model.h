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

// cut selection strategies
#define MOST_VIOLATED_CUT 1
#define ORTHOGONAL_CUTS 2
#define ALL_CUTS 3

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

bool inline sort_pairs_by_snd_val(pair<long,double> a, pair<long,double> b)
{
    return ( a.second < b.second ); 
}

class cut_statistics
{
public:
    map<string,int> pool;               // stores cuts not to add repeated ones

    vector<int> sec_diff_cuts;          // # of different cuts in each separation
    vector<double> sec_infeas_std_dev;  // std.dev. of the amount violated    
};

/// information of a violated subtour elimination constraint
class violated_sec
{
public:
    violated_sec(long edge_count)
    {
        this->edge_count = edge_count;
        coefficients = vector<bool>(edge_count, false);
    }

    virtual ~violated_sec() { }
    
    string toString()
    {
        ostringstream sec_lhs;
        sec_lhs.str("");
        for (vector<bool>::iterator it = coefficients.begin(); it != coefficients.end(); ++it)
            sec_lhs << *it;
        return sec_lhs.str();
    }

     void reset()
     {
        S.clear();
        vertex_count = 0;
        infeasibility = 0.;
        coefficients = vector<bool>(edge_count, false);
     }

    long edge_count;            // instance parameter
    vector<long> S;             // index of edge vars in cut
    long vertex_count;          // number of vertices spanned by S
    double infeasibility;       // amount by which the sec is violated
    vector<bool> coefficients;  // hyperplane coefficients
};

class StableSpanningTreeModel: public KStabModel
{
public:
    StableSpanningTreeModel(IO*);
    virtual ~StableSpanningTreeModel();

    int solve(bool);             // full IP (with SEC in the ORIGINAL GRAPH)

    bool solve_lp_relax(bool);   // corresponding LP relaxation
    double lp_bound;
    double lp_runtime;
    long lp_passes;

private:
    friend class violated_sec;   // encapsulates a class of violated cuts

    cut_statistics *stats;
};

#endif
