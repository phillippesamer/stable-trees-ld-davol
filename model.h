#ifndef _MODEL_H_
#define _MODEL_H_

#include <iostream>
#include <sstream>
#include <cmath>

#include "gurobi_c++.h"

#include "io.h"

#define EPSILON_TOL 0.00001

class Model
{
public:
    Model(IO*);
    ~Model();
    
    int solve();
    int solve_lp_relax();
    bool check_half_integer_solution();

private:
    IO *instance;

    GRBEnv *env;
    GRBModel *model;
    GRBVar *x;

    long card_k;

    void create_objective();
    void create_variables();
    void create_constraints();

    // model for a conflict graph (e.g. from an MSTCC input instance)
    void create_objective_cc();
    void create_variables_cc();
    void create_constraints_cc();
};

#endif
