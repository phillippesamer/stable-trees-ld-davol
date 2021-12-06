#ifndef _MSTCCMODEL_H_
#define _MSTCCMODEL_H_

#include <iostream>
#include <sstream>
#include <sys/time.h>

#include "gurobi_c++.h"

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
    ~StableSpanningTreeModel();
    
    bool solve_lp_relax(bool); // full model (with SEC in the ORIGINAL GRAPH)
    double lp_bound;
    double lp_runtime;

private:
    bool separate_SEC();
};

#endif
