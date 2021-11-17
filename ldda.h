#ifndef _LDDA_H_
#define _LDDA_H_

#include "model.h"
#include <algorithm>

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

    bool dual_ascent(bool);
    vector< vector<long> > multipliers_log;  // lambda[iteration][var_index]
    vector<long> bound_log;
    vector< pair<long, vector<bool> > > solution_pool; // feasible points found

    vector< pair<long,bool> > fixed_vars;   // vars fixed during execution
    IO* flush_fixes_to_instance();

private:
    IO *instance;
    Model *model;

    long edge_deletion_bound();
    long edge_contraction_bound();
    long vertex_deletion_bound();
    long vertex_fix_bound();

    void print_edge_weights();

    vector<long> original_weights;
};

#endif
