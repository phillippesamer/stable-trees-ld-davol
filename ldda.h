#ifndef _LDDA_H_
#define _LDDA_H_

#include "model.h"

/***
 * \file ldda.h
 * 
 * Module for the operations within Lagrangean Decomposition based dual ascent.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 04.11.2021
 */
class LDDA
{
public:
    LDDA(IO*, Model*);
    LDDA(IO*, Model*, vector<double>);
    ~LDDA();

    double dual_ascent();

private:
    vector< vector<double> > multipliers_log;  // lambda[iteration][var_index]

    IO *instance;
    Model *model;

    double edge_deletion_bound();
    double edge_contraction_bound();
    double vertex_deletion_bound();
    double vertex_fix_bound();

    void print_edge_weights();
};

#endif
