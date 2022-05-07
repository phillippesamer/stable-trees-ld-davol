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

#include "io.h"
#include "kstab_model.h"

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
};

#endif
