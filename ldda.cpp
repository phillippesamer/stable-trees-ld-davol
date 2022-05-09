#include "ldda.h"

// algorithm setup switches
bool USE_LDDA_TIME_LIMIT = true;
double LDDA_TIME_LIMIT = 1800;
bool SUBPROBLEM_TIMES_WITHOUT_PROBING = true;  // log information
double TIME_LOG_INTERVAL = 300;

LDDA::LDDA(IO *instance, KStabModel *model)
{
    this->instance = instance;
    this->model = model;
    this->bound_log = vector<double>();
    this->solution_pool = vector< pair<double, vector<bool> > >();
    this->full_log = stringstream();
    this->problem_solved = false;

    this->runtime = -1;
    this->ldda_clock_start = (struct timeval *) malloc(sizeof(struct timeval));
    this->ldda_clock_partial = (struct timeval *) malloc(sizeof(struct timeval));
    this->ldda_clock_stop = (struct timeval *) malloc(sizeof(struct timeval));

    original_weights = vector<double>(instance->graph->w);

    // initialize multipliers at half the original weights
    multipliers_log = vector< vector<double> >();
    multipliers_log.push_back( vector<double>(instance->graph->num_edges, 0) );
    for (long idx=0; idx < instance->graph->num_edges; ++idx)
        multipliers_log[0][idx] = std::round(original_weights[idx] / 2);
}

LDDA::LDDA(IO *instance, KStabModel *model, vector<double> initial_multipliers)
{
    this->instance = instance;
    this->model = model;
    this->bound_log = vector<double>();
    this->solution_pool = vector< pair<double, vector<bool> > >();
    this->full_log = stringstream();
    this->problem_solved = false;

    this->runtime = -1;
    this->ldda_clock_start = (struct timeval *) malloc(sizeof(struct timeval));
    this->ldda_clock_partial = (struct timeval *) malloc(sizeof(struct timeval));
    this->ldda_clock_stop = (struct timeval *) malloc(sizeof(struct timeval));

    original_weights = vector<double>(instance->graph->w);

    // use given values as the initial multipliers
    multipliers_log = vector< vector<double> >();
    multipliers_log.push_back(initial_multipliers);
}

LDDA::~LDDA()
{
    this->multipliers_log.clear();
    this->bound_log.clear();
    this->original_weights.clear();
    this->solution_pool.clear();

    free(ldda_clock_start);
    free(ldda_clock_partial);
    free(ldda_clock_stop);
}

bool LDDA::dual_ascent(bool steepest_ascent)
{
    /***
     * Main loop
     * The objective function in the Lagrangean Decomposition formulation is
     * z(\lambda) = min {(w - \lambda) x} + min {\lambda y}, with x taken in
     * the set of spanning trees in G and y taken in the set of stable sets of
     * cardinality |V|-1 in the conflict graph \hat{G}.
     * Returns true if the execution stopped from a planned outcome, and false
     * if some step failed.
     */

    // log format
    cout << "LDDA: Lagrangean Decomposition bound approximation by Dual Ascent" << endl;

    cout << setw(63) << "";
    cout << setw(9) << "mst  ";
    cout << setw(9) << "mst  ";
    cout << setw(11) << "kstab ";
    cout << setw(11) << "kstab ";
    cout << setw(35) << "" << endl;

    cout << setw(9) << "feasible?";
    cout << setw(9) << "iter";
    cout << setw(9) << "bound";
    cout << setw(12) << "#attempts";
    cout << setw(12) << "direction";
    cout << setw(12) << "adjustment";
    cout << setw(9) << "weight";
    cout << setw(9) << "runtime";
    cout << setw(11) << "weight";
    cout << setw(11) << "runtime";
    cout << setw(13) << "iter (s)";
    cout << setw(7) << "obs";
    cout << endl;

    // procedure timer
    this->start_timer();
    double previous_time_stamp = 0;

    // 0. INITIALIZE SUBPROBLEMS' WEIGHTS WITH LAGRANGEAN MULTIPLIERS

    // args: source1 begin, source1 end, source2 begin, dest begin, operator
    transform( instance->graph->w.begin(),
               instance->graph->w.end(),
               multipliers_log[0].begin(),
               instance->graph->w.begin(),
               minus<double>() );

    instance->graph->update_all_weights(instance->graph->w);
    model->update_all_weights(multipliers_log[0]);

    // iterate while some direction admits a positive step size
    iter = 0;
    problem_solved = false;
    time_limit_exceeded = false;

    do
    {
        ++iter;
        iter_update = false;

        // information for screen log
        bool feasible_mst = false;
        stringstream feasible_mst_msg;

        bool feasible_kstab = false;
        stringstream feasible_kstab_msg;

        // total time spent on solving the (mst/kstab) subproblems and probing 
        double total_mst_time = 0;
        double total_kstab_time = 0;

        // 1. SOLVE MST(G, w-lambda^r)

        if (instance->graph->mst() == false)
        {
            cout << endl << "Graph::mst() failed - graph is not connected"
                 << endl << "infeasible problem instance" << endl;

            this->problem_solved = true;
            this->runtime = total_time();
            return false;
        }

        total_mst_time += instance->graph->mst_runtime;

        // 2. SOLVE KSTAB(\hat{G}, lambda^r)

        if (model->solve(false) <= 0 || model->solution_status != AT_OPTIMUM)
        {
            cout << endl << "KStabModel::solve() failed" 
                 << endl << "infeasible problem instance" << endl;

            this->problem_solved = true;
            this->runtime = total_time();
            return false;
        }

        total_kstab_time += model->runtime();

        // now we can determine the lagrangean bound corresponding to the first multipliers
        if (iter == 1)
            bound_log.push_back(instance->graph->mst_weight + model->solution_weight);

        // 3. COLLECT SET OF INDICES OF VARS WHERE THE SOLUTIONS DON'T MATCH

        vector<long> mismatch;
        for (long idx=0; idx < instance->graph->num_edges; ++idx)
            if (instance->graph->mst_vector[idx] != model->solution_vector[idx])
                mismatch.push_back(idx);

        // 4. CHECK PRIMAL FEASIBILITY OF SOLUTIONS TO SUBPROBLEMS

        // 4.1 THE SET IS EMPTY: THE SOLUTIONS MATCH AND THE PROBLEM IS SOLVED
        if (mismatch.empty())
        {
            cout << "the solutions to both subproblems are the same" << endl;
            cout << "problem solved to optimality" << endl;

            problem_solved = true;

            this->runtime = total_time();
            return true;
        }

        // 4.2 IF THE MST SOLUTION IS STABLE, IT IS AN INTEGER FEASIBLE POINT
        if ( instance->test_stability(instance->graph->mst_vector) )
        {
            // store solution and its cost (wrt original weights)
            double true_cost = 0;
            double cost_in_other_subproblem = 0;
            for (long idx=0; idx < instance->graph->num_edges; ++idx)
            {
                // if x(idx) = 1
                if (instance->graph->mst_vector[idx])
                {
                    true_cost += original_weights[idx];
                    cost_in_other_subproblem += multipliers_log.back()[idx];
                }
            }

            this->solution_pool.push_back( make_pair(true_cost, instance->graph->mst_vector) );

            feasible_mst = true;
            feasible_mst_msg = stringstream();
            feasible_mst_msg << "mst primal feasible (weight " << true_cost << ")";

            // if it has the same lambda^r cost as the kstab, than it is optimal
            if (cost_in_other_subproblem == model->solution_weight)
            {
                cout << "the mst solution is primal feasible and dual optimal ; "
                     << "problem solved to optimality" << endl;

                problem_solved = true;

                this->runtime = total_time();
                return true;
            }
        }

        // 4.3 IF THE KSTAB SOLUTION IS ACYCLIC, IT IS AN INTEGER FEASIBLE POINT
        if ( instance->test_acyclic_kstab(model->solution_vector) )
        {
            // store solution and its cost (wrt original weights)
            double true_cost = 0;
            double cost_in_other_subproblem = 0;
            for (long idx=0; idx < instance->graph->num_edges; ++idx)
            {
                // if y(idx) = 1
                if (model->solution_vector[idx])
                {
                    true_cost += original_weights[idx];
                    cost_in_other_subproblem += ( original_weights[idx]
                                                - multipliers_log.back()[idx] );
                }
            }

            this->solution_pool.push_back( make_pair(true_cost, model->solution_vector) );

            feasible_kstab = true;
            feasible_kstab_msg = stringstream();
            feasible_kstab_msg << "kstab primal feasible (weight " << true_cost << ")";

            // if it has the same (w-lambda^r) cost as the mst, than it is optimal
            if (cost_in_other_subproblem == instance->graph->mst_weight)
            {
                cout << "the kstab solution is primal feasible and dual optimal ; "
                     << "problem solved to optimality" << endl;

                problem_solved = true;

                this->runtime = total_time();
                return true;
            }
        }

        /* 5. UNLESS WE SOLVED THE PROBLEM, CHOOSE ELEMENT FOR DUAL ASCENT
         *
         * STEEPEST ASCENT: evaluate each ascent direction, and follow the one
         * yielding the greatest bound.
         * 
         * FIRST ASCENT: follow the first maximal ascent direction available.
         * To impose an order on the set of mismatching variables (and seeking
         * not to favour some of them throughout execution), we start indexing
         * from the iteration number (and mod out by the cardinality of the
         * mismatch set).
         */

        long chosen_direction = -1;
        double chosen_adjustment = 0;
        double chosen_bound_improvement = -1;

        unsigned long attempt = 0;
        long attempting_idx = (iter-1) % mismatch.size();  // starting direction

        bool done_with_first_ascent = false;

        while ( !done_with_first_ascent     && 
                attempt < mismatch.size()   &&
                !time_limit_exceeded        )
        {
            long current_direction = mismatch[attempting_idx];

            // case 1: x^r_e = 1 and y^r_e = 0 (inoc2022, thm 4.3)
            if (instance->graph->mst_vector[current_direction])
            {
                #ifdef DEBUG
                    /*
                    cout << "x^" << iter << "_" << current_direction 
                         << " = 1 > 0 = y^"  << iter << "_"
                         << current_direction << endl;
                    */

                    if (model->solution_vector[current_direction])
                    {
                        cout << "UNEXPECTED ERROR: ";
                        cout << "x^" << iter << "_" << current_direction
                             << "  =  1  =  y^"  << iter << "_"
                             << current_direction << endl;
                        this->runtime = total_time();
                        return false;
                    }
                #endif

                // PROBING MST WITHOUT e TO DETERMINE \del^r_e in Thm 4.3 (INOC)
                pair<bool,double> probing_mst = instance->graph->mst_probing_var(current_direction, false);
                total_mst_time += instance->graph->probe_runtime;

                // PROBING KSTAB FORCING e TO DETERMINE \delta^r_e in Thm 4.3 (INOC)
                pair<ModelStatus,double> probing_kstab = model->probe_var(current_direction, true);
                total_kstab_time += model->runtime();

                if (probing_kstab.first == STATUS_UNKNOWN)
                {
                    cout << "probing var y[" << current_direction << "] = 1 failed" 
                         << " (runtime: " << model->runtime() << " s)" << endl;

                    this->runtime = total_time();
                    return false;
                }

                if (probing_kstab.first == IS_INFEASIBLE && probing_mst.first == false)
                {
                    // a mismatching variable cannot agree in a feasible solution => infeasible instance
                    cout << endl << "mismatching variable (index " << current_direction
                         << ") cannot agree on both subproblems" << endl
                         << "infeasible problem instance" << endl;

                    this->problem_solved = true;
                    this->runtime = total_time();
                    return false;
                }

                /***
                 * AT LEAST ONE OF THE PROBES IS FEASIBLE => PROCEED TO COMPUTE THE BOUNDS
                 * In this case, we set the infeasible probe value to infinity, so that
                 * an ascent in this direction is possible iff the other probe is
                 * feasible and the bound increase (del/delta) is positive
                 */

                double del = (probing_mst.first == true) ?
                    probing_mst.second - instance->graph->mst_weight
                    : numeric_limits<double>::max();

                double delta = (probing_kstab.first == AT_OPTIMUM) ?
                    probing_kstab.second - model->solution_weight
                    : numeric_limits<double>::max();

                #ifdef DEBUG
                    if (min(del,delta) < 0)
                    {
                        cout << "UNEXPECTED ERROR: min( delta=" << delta << ", del=" << del <<" ) < 0" << endl;
                        this->runtime = total_time();
                        return false;
                    }
                    /*
                    else
                        cout << "delta = " << delta << ", del = " << del << endl;
                    */
                #endif

                if ( min(del,delta) > 0 &&
                     min(del,delta) > chosen_bound_improvement )
                {
                    iter_update = true;

                    chosen_direction = current_direction;
                    chosen_adjustment = (-1)*min(del,delta);
                    chosen_bound_improvement = min(del,delta);

                    if (!steepest_ascent)
                        done_with_first_ascent = true;
                }

            }

            // case 2: x^r_e = 0 and y^r_e = 1 (inoc2022, thm 4.2)
            else
            {
                #ifdef DEBUG
                    /*
                    cout << "x^" << iter << "_" << current_direction 
                         << " = 0 < 1 = y^"  << iter << "_"
                         << current_direction << endl;
                    */

                    if (!model->solution_vector[current_direction])
                    {
                        cout << "UNEXPECTED ERROR: ";
                        cout << "x^" << iter << "_" << current_direction
                             << "  =  0  =  y^"  << iter << "_"
                             << current_direction << endl;

                        this->runtime = total_time();
                        return false;
                    }
                #endif

                // PROBING MST FORCING e TO DETERMINE \del^r_e in Thm 4.2 (INOC)
                pair<bool,double> probing_mst = instance->graph->mst_probing_var(current_direction, true);
                total_mst_time += instance->graph->probe_runtime;

                // contracting an edge does not make a connected graph disconnected
                if (probing_mst.first == false)
                {
                    cout << "UNEXPECTED ERROR: ";
                    cout << "probing x^" << iter << "_" << current_direction
                         << "  =  1 (edge contraction) gave disconnected graph"
                         << endl;

                    this->runtime = total_time();
                    return false;
                }

                // PROBING KSTAB WITHOUT e TO DETERMINE \delta^r_e in Thm 4.2 (INOC)
                pair<ModelStatus,double> probing_kstab = model->probe_var(current_direction, false);
                total_kstab_time += model->runtime();

                // the call above took care of the case where the probe is infeasible
                if (probing_kstab.first == STATUS_UNKNOWN)
                {
                    cout << "probing var y[" << current_direction << "] = 0 failed" 
                         << " (runtime: " << model->runtime() << " s)" << endl;

                    this->runtime = total_time();
                    return false;
                }

                /***
                 * AT LEAST ONE OF THE PROBES IS FEASIBLE => PROCEED TO COMPUTE THE BOUNDS
                 * In this case, we set the infeasible probe value to infinity, so that
                 * an ascent in this direction is possible iff the other probe is
                 * feasible and the bound increase (del/delta) is positive
                 */

                double del = probing_mst.second - instance->graph->mst_weight;

                double delta = (probing_kstab.first == AT_OPTIMUM) ?
                    probing_kstab.second - model->solution_weight
                    : numeric_limits<double>::max();

                #ifdef DEBUG
                    if (min(del,delta) < 0)
                    {
                        cout << "UNEXPECTED ERROR: min( delta=" << delta << ", del=" << del <<" ) < 0" << endl;
                        this->runtime = total_time();
                        return false;
                    }
                    /*
                    else
                        cout << "delta = " << delta << ", del = " << del << endl;
                    */
                #endif

                if ( min(del,delta) > 0 &&
                     min(del,delta) > chosen_bound_improvement )
                {
                    iter_update = true;

                    chosen_direction = current_direction;
                    chosen_adjustment = min(del,delta);
                    chosen_bound_improvement = min(del,delta);

                    if (!steepest_ascent)
                        done_with_first_ascent = true;
                }

            }

            // done evaluating this direction (i.e. mismatching variable)

            attempting_idx = (attempting_idx+1) % mismatch.size();
            ++attempt;

            if (USE_LDDA_TIME_LIMIT && total_time() > LDDA_TIME_LIMIT)
                time_limit_exceeded = true;
        }

        // 6. UPDATE MULTIPLIER \lambda^{r+1}_e, COPY THE OTHERS

        if (chosen_direction >= 0)
        {
            vector<double> next_multipliers = vector<double>( multipliers_log.back() );
            next_multipliers[chosen_direction] = next_multipliers[chosen_direction]
                                                 + chosen_adjustment;
            multipliers_log.push_back(next_multipliers);

            double new_bound = bound_log.back() + chosen_bound_improvement;
            bound_log.push_back(new_bound);

            // update objective in mst subproblem: w_e - lambda_e
            instance->graph->update_single_weight(chosen_direction,
                original_weights[chosen_direction] - next_multipliers[chosen_direction]);

            // update objective in kstab subproblem: lambda_e
            model->update_single_weight(chosen_direction,
                next_multipliers[chosen_direction]);
        }

        // iteration time
        double iter_time = partial_time();

        // 7. SCREEN LOG
        stringstream logline;
        logline.precision(4);

        if (feasible_mst)
            logline << " * ";
        else
            logline << " - ";

        if (feasible_kstab)
            logline << " * ";
        else
            logline << " - ";

        logline << " " << setw(9) << iter;
        logline << " " << setw(9) << bound_log.at(iter-1);
        logline << " " << setw(9) << attempt << "/" << mismatch.size();
        if (chosen_direction >= 0)
        {
            logline << " " << setw(9) << chosen_direction;
            logline << " " << setw(9) << chosen_adjustment;
        }
        else
        {
            logline << " " << setw(9) << "-";
            logline << " " << setw(9) << "-";
        }
        logline << " " << setw(11) << instance->graph->mst_weight;

        if (SUBPROBLEM_TIMES_WITHOUT_PROBING)
            logline << " " << fixed << setw(10) << instance->graph->mst_runtime;
        else
            logline << " " << fixed << setw(10) << total_mst_time;

        logline << " " << setw(9) << model->solution_weight;

        if (SUBPROBLEM_TIMES_WITHOUT_PROBING)
            logline << " " << fixed << setw(11) << model->solution_runtime;
        else
            logline << " " << fixed << setw(11) << total_kstab_time;

        logline << " " << fixed << setw(12) << iter_time;

        //logline << right << setw(0);
        //logline.unsetf(ios_base::floatfield);

        if (feasible_mst || feasible_kstab)
            logline << setw(6) << "";

        if (feasible_mst)
            logline << feasible_mst_msg.str();

        if (feasible_mst && feasible_kstab)
            logline << "; ";

        if (feasible_kstab)
            logline << feasible_kstab_msg.str();

        if (time_limit_exceeded)
            logline << " time limit exceeded";
        else
        {
            // log total time at set intervals
            double current_elapsed_time = total_time();
            if (current_elapsed_time >= previous_time_stamp + TIME_LOG_INTERVAL)
            {
                logline << " ldda total time: " << current_elapsed_time;
                previous_time_stamp = current_elapsed_time;
            }
        }

        cout << logline.str() << endl;      // screen log
        full_log << logline.str() << endl;  // full log, dumped when create_log() is called
    }
    while (iter_update && !time_limit_exceeded);

    this->runtime = total_time();

    return true;
}

stringstream LDDA::create_log()
{
    /// prepare log of the LDDA execution as a stringstream object
    stringstream log;

    log << "LDDA: Lagrangean Decomposition bound approximation by Dual Ascent" << endl;

    log << setw(63) << "";
    log << setw(9) << "mst  ";
    log << setw(9) << "mst  ";
    log << setw(11) << "kstab ";
    log << setw(11) << "kstab ";
    log << setw(35) << "" << endl;

    log << setw(9) << "feasible?";
    log << setw(9) << "iter";
    log << setw(9) << "bound";
    log << setw(12) << "#attempts";
    log << setw(12) << "direction";
    log << setw(12) << "adjustment";
    log << setw(9) << "weight";
    log << setw(9) << "runtime";
    log << setw(11) << "weight";
    log << setw(11) << "runtime";
    log << setw(13) << "iter (s)";
    log << setw(7) << "obs";
    log << endl;

    log << full_log.str() << endl;

    log << "LDDA multipliers log:" << endl;
    for (unsigned i=0; i < multipliers_log.size(); ++i)
    {
        log << "#"<< setw(3) << left;
        log << i;
        log << setw(0) << right;
        log << ": ";
        for (unsigned idx=0; idx < multipliers_log[i].size(); ++idx)
        {
            //if (idx+1 % 30 == 0) log << endl;
            log << multipliers_log[i][idx] << " ";
        }
        log << endl;
    }

    return log;
}

void LDDA::start_timer()
{
    gettimeofday(ldda_clock_start, 0);
    gettimeofday(ldda_clock_partial, 0);  // initialization for the 1st call to partial_time()
}

double LDDA::partial_time()
{
    // returns the time elapsed since last call
    struct timeval* tmp = (struct timeval *) malloc(sizeof(struct timeval));
    gettimeofday(tmp, 0);

    unsigned long clock_time = 1.e6 * (tmp->tv_sec - ldda_clock_partial->tv_sec) +
                                      (tmp->tv_usec - ldda_clock_partial->tv_usec);

    free(tmp);
    gettimeofday(ldda_clock_partial, 0);

    return (double) clock_time / (double)1.e6 ;
}

double LDDA::total_time()
{
    gettimeofday(ldda_clock_stop, 0);

    unsigned long clock_time = 1.e6 * (ldda_clock_stop->tv_sec - ldda_clock_start->tv_sec) +
                                      (ldda_clock_stop->tv_usec - ldda_clock_start->tv_usec);

    return (double) clock_time / (double)1.e6 ;
}
