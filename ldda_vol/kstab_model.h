#ifndef _KSTABMODEL_H_
#define _KSTABMODEL_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>

#include "gurobi_c++.h"

#include "io.h"

#define EPSILON_TOL 0.00000001
#define ZERO_TOL 0.0001

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
    double solution_weight;
    double solution_dualbound;
    vector<bool> solution_vector;
    ModelStatus solution_status;
    double solution_runtime;

    double runtime();
    void set_time_limit(double);

    void update_single_weight(long,double);
    void update_all_weights(vector<double>);

    pair<ModelStatus,double> probe_var(long, bool);

    vector<long> fix_var(long, bool);

protected:
    friend class LDDA;
    friend class LDDAVolume;

    IO *instance;
    long fixed_cardinality;

    GRBEnv *env;
    GRBModel *model;
    GRBVar *x;

    // model for kstabs IN THE CONFLICT GRAPH (in an MSTCC input instance)
    void create_objective();
    void create_variables();
    void create_constraints();

    int save_optimization_status();

    // structures for enumerating maximal cliques in the conflict graph
    long clique_counter;
    map<long,long> clique_sizes;
    vector<vector<bool> > cliques_bit_adj;
    void all_maximal_cliques(vector<bool> &, vector<bool> &, GRBLinExpr &);
};

#endif
