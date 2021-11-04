#ifndef _MODEL_H_
#define _MODEL_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>

#include "gurobi_c++.h"

#include "io.h"

#define EPSILON_TOL 0.00001

enum ProbeStatus {PROBE_OPTIMAL, PROBE_INFEASIBLE, PROBE_UNKNOWN};

/***
 * \file model.h
 * 
 * Module for the integer programming model to find fixed cardinality stable
 * sets of minimum weight, using the Gurobi solver API.
 * 
 * NB! This implementation finds fixed cardinality stable sets in the conflict
 * graph of the corresponding MSTCC instance.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 02.11.2021
 */
class Model
{
public:
    Model(IO*);
    ~Model();
    
    int solve();
    double runtime();

    pair<ProbeStatus,double> probe_var_at_zero(long);
    pair<ProbeStatus,double> probe_var_at_one(long);

private:
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
