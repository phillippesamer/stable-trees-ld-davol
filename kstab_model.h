#ifndef _KSTABMODEL_H_
#define _KSTABMODEL_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>

#include "gurobi_c++.h"

#include "io.h"

#define EPSILON_TOL 0.00000001

enum ModelStatus {AT_OPTIMUM, IS_INFEASIBLE, STATUS_UNKNOWN};

/***
 * \file kstab_model.h
 * 
 * Module for the integer programming model to find fixed cardinality stable
 * sets of minimum weight, using the Gurobi solver API.
 * 
 * NB! This implementation finds fixed cardinality stable sets in the conflict
 * graph of the corresponding MSTCC instance.
 * 
 * LDDA class is declared a friend to avoid cumbersome get/set calls.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 02.11.2021
 */
class KStabModel
{
public:
    KStabModel(IO*);
    virtual ~KStabModel();
    
    int solve(bool);
    long solution_weight;
    vector<bool> solution_vector;
    ModelStatus solution_status;
    double solution_runtime;

    double runtime();

    void update_single_weight(long,long);
    void update_all_weights(vector<long>);

    pair<ModelStatus,long> probe_var(long, bool);

    vector<long> fix_var(long, bool);

protected:
    friend class LDDA;

    IO *instance;
    long fixed_cardinality;

    GRBEnv *env;
    GRBModel *model;
    GRBVar *x;

    // model for kstabs IN THE CONFLICT GRAPH (in an MSTCC input instance)
    void create_objective();
    void create_variables();
    void create_constraints();
};

#endif
