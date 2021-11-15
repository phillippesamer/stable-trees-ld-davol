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
    LDDA(IO*, Model*, vector<long>);
    ~LDDA();

    long dual_ascent(bool);

private:
    vector< vector<long> > multipliers_log;  // lambda[iteration][var_index]

    IO *instance;
    Model *model;

    long edge_deletion_bound();
    long edge_contraction_bound();
    long vertex_deletion_bound();
    long vertex_fix_bound();

    void print_edge_weights();

    vector< pair<long,bool> > fixed_vars;   // vars fixed during execution
};

#endif
