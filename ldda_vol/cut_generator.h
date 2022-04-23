#ifndef _MSTCC_CUT_GEN_H_
#define _MSTCC_CUT_GEN_H_

#include <iostream>
#include <sstream>
#include <iomanip>

#include "gurobi_c++.h"

#include "io.h"
#include "kstab_model.h"
#include "mstcc_model.h"

/***
 * \file cut_generator.h
 * 
 * Module for the Gurobi callback class, implementing the separation procedure
 * for subtour elimination constraints (SEC).
 * 
 * Extends (and implements) the abstract base class in Gurobi.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 29.03.2022
 */

class CutGenerator: public GRBCallback
{
public:
    CutGenerator(GRBVar*, IO*, cut_statistics*);
    virtual ~CutGenerator();

protected:
    void callback();
    void separate_lpr();

private:
    friend class KStabModel;
    friend class StableSpanningTreeModel;

    // input and model data
    IO *instance;
    GRBModel *model;
    cut_statistics *stats;

    GRBVar* x_vars;
    double *x_val;
    long num_vars;

    long sec_counter;
    void separate_sec(int);

    bool separate_SEC_integer(vector<GRBLinExpr> &, vector<long> &);
    bool separate_SEC_folklore(vector<GRBLinExpr> &, vector<long> &);
    bool separate_SEC_classical(vector<GRBLinExpr> &, vector<long> &);
    void dfs_checking_acyclic(long, long, vector<bool> &, long &, stack<long> &, vector< list<long> > &, bool &);
};

#endif
