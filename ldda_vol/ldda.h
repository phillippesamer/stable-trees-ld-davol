#ifndef _LDDA_H_
#define _LDDA_H_

#include "kstab_model.h"
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <sys/time.h>

/***
 * \file ldda.h
 * 
 * Module for computing (approximately) the Lagrangean Decomposition bound
 * with a dual ascent algorithm.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 04.11.2021
 */
class LDDA
{
public:
    LDDA(IO*, KStabModel*);
    LDDA(IO*, KStabModel*, vector<double>);
    virtual ~LDDA();

    bool dual_ascent(bool);
    bool problem_solved;
    bool iter_update;
    bool time_limit_exceeded;
    long iter;
    double runtime;

    vector< vector<double> > multipliers_log;  // lambda[iteration][var_index]
    vector<double> bound_log;
    vector< pair<double, vector<bool> > > solution_pool; // (weight, feasible point)

    vector< pair<long,bool> > fixed_vars;   // vars fixed during execution
    IO* flush_fixes_to_instance();

    stringstream create_log();

protected:
    IO *instance;
    KStabModel *model;

    pair<bool,double> edge_deletion_bound(long);
    pair<bool,double> edge_contraction_bound(long);
    pair<ModelStatus,double> vertex_deletion_bound(long);
    pair<ModelStatus,double> vertex_fix_bound(long);
    void fix_element_at_one_in_graph_and_model(long);

    vector<double> original_weights;

    double contracted_edges_weight;
    vector<long> contracted_edges;
    vector<bool> contracted_edges_mask;
    vector<bool> removed_edges_mask;

    stringstream full_log;

    void start_timer();
    double partial_time();
    double total_time();
    struct timeval* ldda_clock_start;
    struct timeval* ldda_clock_partial;
    struct timeval* ldda_clock_stop;
};

#endif
