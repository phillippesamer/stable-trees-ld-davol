#include "ldda.h"

LDDA::LDDA(IO *instance, Model *model)
{
    this->instance = instance;
    this->model = model;
    this->fixed_vars = vector< pair<long,bool> >();
    this->bound_log = vector<long>();
    this->solution_pool = vector< pair<long, vector<bool> > >();

    this->contracted_edges_weight = 0;
    this->contracted_edges = vector<long>();
    this->contracted_edges_mask = vector<bool>(instance->graph->num_edges, false);

    // backup original weights to restore later
    original_weights = vector<long>(instance->graph->w);

    // initialize multipliers at zero
    multipliers_log = vector< vector<long> >();
    multipliers_log.push_back( vector<long>(instance->graph->num_edges, 0) );
}

LDDA::LDDA(IO *instance, Model *model, vector<long> initial_multipliers)
{
    this->instance = instance;
    this->model = model;
    this->fixed_vars = vector< pair<long,bool> >();
    this->bound_log = vector<long>();
    this->solution_pool = vector< pair<long, vector<bool> > >();

    this->contracted_edges_weight = 0;
    this->contracted_edges = vector<long>();
    this->contracted_edges_mask = vector<bool>(instance->graph->num_edges, false);

    // backup original weights to restore later
    original_weights = vector<long>(instance->graph->w);

    // use given values as the initial multipliers
    multipliers_log = vector< vector<long> >();
    multipliers_log.push_back(initial_multipliers);
}

LDDA::~LDDA()
{
    this->multipliers_log.clear();
    this->fixed_vars.clear();
    this->bound_log.clear();
    this->original_weights.clear();
    this->solution_pool.clear();
    this->contracted_edges.clear();
    this->contracted_edges_mask.clear();
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

    // 0. INITIALIZE WEIGHTS WITH LAGRANGEAN MULTIPLIERS

    // args: source1 begin, source1 end, source2 begin, dest begin, operator
    transform( instance->graph->w.begin(),
               instance->graph->w.end(),
               multipliers_log[0].begin(),
               instance->graph->w.begin(),
               minus<long>() );

    #ifdef DEBUG
        cout << "initial mst objective vector: ";
        print_edge_weights();
    #endif

    instance->graph->update_all_weights(instance->graph->w);

    model->update_all_weights(multipliers_log[0]);

    iter = 1;
    problem_solved = false;

    do
    {
        iter_update = false;

        // 1. SOLVE MST(G, w-lambda^r)
        if (instance->graph->mst() == false)
        {
            cout << "Graph::mst() failed (graph is not connected)" << endl;
            cout << "infeasible problem instance" << endl;
            return false;
        }

        /* For separation of concerns purposes, edges contracted by LDDA (in
         * the LEMON data structure only!) are missing in later spanning trees,
         * so we include them here
         */
        instance->graph->mst_weight += contracted_edges_weight;
        for ( vector<long>::iterator it = contracted_edges.begin();
              it != contracted_edges.end(); ++it )
        {
            instance->graph->mst_vector[*it] = true;
        }

        // 2. SOLVE KSTAB(\hat{G}, lambda^r)
        if (model->solve(true) <= 0 || model->solution_status != AT_OPTIMUM)
        {
            cout << "Model::solve() failed" << endl;
            cout << "infeasible problem instance" << endl;
            return false;
        }

        // 3. COLLECT SET OF INDICES OF VARS WHERE THE SOLUTIONS DON'T MATCH
        vector<long> mismatch;
        for (long idx=0; idx < instance->graph->num_edges; ++idx)
        {
            if (instance->graph->mst_vector[idx] != model->solution_vector[idx])
            {
                mismatch.push_back(idx);

                #ifdef DEBUG
                    if (contracted_edges_mask[idx])
                    {
                        cout << "UNEXPECTED ERROR: solutions differ in an "
                             << "index corresponding to a contracted edge"
                             << endl;
                        return false;
                    }
                #endif
            }
        }

        // 4. CHECK PRIMAL FEASIBILITY OF SOLUTIONS TO SUBPROBLEMS

        // 4.1 THE SET IS EMPTY: THE SOLUTIONS MATCH AND THE PROBLEM IS SOLVED
        if (mismatch.empty())
        {
            cout << "the solutions to both subproblems are the same" << endl;
            cout << "problem solved to optimality" << endl;

            problem_solved = true;

            return true;
        }

        // 4.2 IF THE MST SOLUTION IS STABLE, IT IS AN INTEGER FEASIBLE POINT
        if ( instance->test_stability(instance->graph->mst_vector) )
        {
            // store solution and its cost (wrt original weights)
            long true_cost = 0;
            long cost_in_other_subproblem = 0;
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

            cout << "the mst solution is primal feasible" << endl;
            cout << "integer feasible point of weight " << true_cost << endl;

            // if it has the same lambda^r cost as the kstab, than it is optimal
            if (cost_in_other_subproblem == model->solution_weight)
            {
                cout << "the mst solution is primal feasible and dual optimal" << endl;
                cout << "problem solved to optimality" << endl;

                problem_solved = true;

                return true;
            }
        }

        // 4.3 IF THE KSTAB SOLUTION IS ACYCLIC, IT IS AN INTEGER FEASIBLE POINT
        if ( instance->test_acyclic_kstab(model->solution_vector) )
        {
            // store solution and its cost (wrt original weights)
            long true_cost = 0;
            long cost_in_other_subproblem = 0;
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

            cout << "the kstab solution is primal feasible" << endl;
            cout << "integer feasible point of weight " << true_cost << endl;

            // if it has the same (w-lambda^r) cost as the mst, than it is optimal
            if (cost_in_other_subproblem == instance->graph->mst_weight)
            {
                cout << "the kstab solution is primal feasible and dual optimal" << endl;
                cout << "problem solved to optimality" << endl;

                problem_solved = true;

                return true;
            }
        }

        // 5. UNLESS WE SOLVED THE PROBLEM, CHOOSE ELEMENT FOR DUAL ASCENT

        long chosen_direction = -1;
        long chosen_adjustment = -1;
        long chosen_bound_improvement = -1;

        if (steepest_ascent)
        {
            /* 5.1 STEEPEST ASCENT
             * Evaluate each ascent direction, and follow the one yielding 
             * the greatest bound
             */

            for (vector<long>::iterator current_direction = mismatch.begin();
                 current_direction != mismatch.end(); ++current_direction)
            {
                // case 1: x^r_e = 1 and y^r_e = 0 (inoc2022, thm 4.3)
                if (instance->graph->mst_vector[*current_direction])
                {
                    #ifdef DEBUG
                        cout << "x^" << iter << "_" << *current_direction 
                             << " = 1 > 0 = y^"  << iter << "_"
                             << *current_direction << endl;

                        if (model->solution_vector[*current_direction])
                        {
                            cout << "UNEXPECTED ERROR: ";
                            cout << "x^" << iter << "_" << *current_direction
                                 << "  =  1  =  y^"  << iter << "_"
                                 << *current_direction << endl;
                            return false;
                        }
                    #endif

                    // PROBING MST WITHOUT e TO DETERMINE \del^r_e
                    pair<bool,long> probing_mst = edge_deletion_bound(*current_direction);

                    // the call above took care of the case where the probe is infeasible
                    if (probing_mst.first == true)
                    {
                        // PROBING KSTAB FORCING e TO DETERMINE \delta^r_e
                        pair<ModelStatus,long> probing_kstab
                            = vertex_fix_bound(*current_direction);

                        // the call above took care of the case where the probe is infeasible
                        if (probing_kstab.first == STATUS_UNKNOWN)
                        {
                            cout << "probing var y[" << *current_direction << "] = 1 failed" 
                            << " (runtime: " << model->runtime() << " s)" << endl;

                            return false;
                        }
                        else if (probing_kstab.first == AT_OPTIMUM)
                        {
                            // probings found feasible solutions => proceed to compute the bounds

                            // NB! contracted edges do not appear in Graph::mst_probing_var()
                            long probing_mst_bound = probing_mst.second + contracted_edges_weight;
                            long probing_kstab_bound = probing_kstab.second;

                            long del = probing_mst_bound - instance->graph->mst_weight;
                            long delta = probing_kstab_bound - model->solution_weight;

                            #ifdef DEBUG
                                if (min(del,delta) < 0)
                                {
                                    cout << "UNEXPECTED ERROR: min( delta=" << delta << ", del=" << del <<" ) < 0" << endl;
                                    return false;
                                }
                            #endif

                            if ( min(del,delta) > 0 &&
                                 min(del,delta) > chosen_bound_improvement )
                            {
                                iter_update = true;

                                chosen_direction = *current_direction;
                                chosen_adjustment = (-1)*min(del,delta);
                                chosen_bound_improvement = min(del,delta);
                            }
                        }
                    }
                }

                // case 2: x^r_e = 0 and y^r_e = 1 (inoc2022, thm 4.2)
                else
                {
                    #ifdef DEBUG
                        cout << "x^" << iter << "_" << *current_direction 
                             << " = 0 < 1 = y^"  << iter << "_"
                             << *current_direction << endl;

                        if (!model->solution_vector[*current_direction])
                        {
                            cout << "UNEXPECTED ERROR: ";
                            cout << "x^" << iter << "_" << *current_direction
                                 << "  =  0  =  y^"  << iter << "_"
                                 << *current_direction << endl;
                            return false;
                        }
                    #endif

                    // PROBING MST FORCING e TO DETERMINE \del^r_e
                    pair<bool,long> probing_mst = edge_contraction_bound(*current_direction);

                    // contracting an edge does not make a connected graph disconnected
                    if (probing_mst.first == false)
                    {
                        cout << "UNEXPECTED ERROR: ";
                        cout << "probing x^" << iter << "_" << *current_direction
                             << "  =  1 (edge contraction) gave disconnected graph"
                             << endl;
                        return false;
                    }

                    // PROBING KSTAB WITHOUT e TO DETERMINE \delta^r_e
                    pair<ModelStatus,long> probing_kstab
                        = vertex_deletion_bound(*current_direction);

                    // the call above took care of the case where the probe is infeasible
                    if (probing_kstab.first == STATUS_UNKNOWN)
                    {
                        cout << "probing var y[" << *current_direction << "] = 0 failed" 
                        << " (runtime: " << model->runtime() << " s)" << endl;

                        return false;
                    }
                    else if (probing_kstab.first == AT_OPTIMUM)
                    {
                        // probings found feasible solutions => proceed to compute the bounds

                        // NB! contracted edges do not appear in Graph::mst_probing_var()
                        long probing_mst_bound = probing_mst.second + contracted_edges_weight;
                        long probing_kstab_bound = probing_kstab.second;

                        long del = probing_mst_bound - instance->graph->mst_weight;
                        long delta = probing_kstab_bound - model->solution_weight;

                        #ifdef DEBUG
                            if (min(del,delta) < 0)
                            {
                                cout << "UNEXPECTED ERROR: min( delta=" << delta << ", del=" << del <<" ) < 0" << endl;
                                return false;
                            }
                        #endif

                        if ( min(del,delta) > 0 &&
                             min(del,delta) > chosen_bound_improvement )
                        {
                            iter_update = true;

                            chosen_direction = *current_direction;
                            chosen_adjustment = min(del,delta);
                            chosen_bound_improvement = min(del,delta);
                        }
                    }
                }

                // done evaluating this direction (= mismatching variable)
            }

        }
        else
        {
            /* 5.2 "FIRST" ASCENT
             * Follow the first maximal ascent direction available. To impose
             * an order on the set of mismatching variables (and seeking not to
             * favour some throughout execution), we use the iteration number
             * mod the cardinality of the mismatch set.
             */

            //start at iter % mismatch.size()
            //do
            // { evaluate step size in this direction } while(not found and not out of directions)
            //save direction and step size

            // switch between x^r_e = 1 and y^r_e = 0 and vice versa
            // determine delta^r_e and del^r_e
            // if probing finds an infeasible problem, fix corresponding variable (save this on fixed_vars), pick new one
        }

        if (chosen_direction < 0 && !iter_update)
        {
            // no direction admits a positive step size - the procedure is over
            cout << "none of the " << mismatch.size() << " directions admits a "
                 << "positive step size"<< endl;
            cout << "dual ascent done" << endl;

            return true;
        }

        // 6. UPDATE MULTIPLIER \lambda^{r+1}_e, COPY THE OTHERS

        // 7. SCREEN LOG
    }
    while (!problem_solved && iter_update);

    if (problem_solved)
    {
        // original problem is solved
    }
    else
    {
        // original problem not solved, but no ascent was possible
    }

    return true;
}

pair<bool,long> LDDA::edge_deletion_bound(long idx)
{
    /// probing mst without given edge to determine \del^r_e in Thm 4.3 (INOC)

    pair<bool,long> probing_mst = instance->graph->mst_probing_var(idx, false);

    // probing var at 0 infeasible => fix var at 1 
    if (probing_mst.first == false)
    {
        cout << "probing var x[" << idx << "] = 0 gives a disconnected graph ";
        cout << "=> fixing element at one throughout" << endl;

        fix_element_at_one_in_graph_and_model(idx);
    }

    return probing_mst;
}

pair<bool,long> LDDA::edge_contraction_bound(long idx)
{
    /// probing mst forcing given edge to determine \del^r_e in Thm 4.2 (INOC)

    pair<bool,long> probing_mst = instance->graph->mst_probing_var(idx, true);

    return probing_mst;
}

pair<bool,long> LDDA::vertex_deletion_bound(long idx)
{
    /// probing kstab without a vertex to determine \delta^r_e in Thm 4.2 (INOC)

    pair<ModelStatus,long> probing_kstab = model->probe_var(idx, false);

    // probing var at 0 infeasible => fix var at 1
    if (probing_kstab.first == IS_INFEASIBLE)
    {
        cout << "probing var y[" << idx << "] = 0 gives an infeasible model "
             << "(runtime: " << model->runtime() << " s)";
        cout << " => fixing element at one throughout" << endl;

        fix_element_at_one_in_graph_and_model(idx);
    }

    return probing_kstab;
}

pair<bool,long> LDDA::vertex_fix_bound(long idx)
{
    /// probing kstab forcing a vertex to determine \delta^r_e in Thm 4.3 (INOC)

    pair<ModelStatus,long> probing_kstab = model->probe_var(idx, true);

    // probing var at 1 infeasible => fix var at 0
    if (probing_kstab.first == IS_INFEASIBLE)
    {
        iter_update = true;
        fixed_vars.push_back( make_pair(idx, false) );

        cout << "probing var y[" << idx << "] = 1 gives an infeasible model "
             << "(runtime: " << model->runtime() << " s)";
        cout << " => fixing element at zero throughout" << endl;

        model->fix_var(idx, false);
        instance->graph->lemon_delete_edge(idx);
    }

    return probing_kstab;
}

void LDDA::fix_element_at_one_in_graph_and_model(long idx)
{
    /***
     * Implement an edge contraction and model variable fix at 1, handling
     * possible implications. Used by both edge_deletion_bound() and
     * vertex_deletion_bound().
     */

    iter_update = true;
    fixed_vars.push_back( make_pair(idx, true) );

    contracted_edges_weight += original_weights[idx];
    contracted_edges.push_back(idx);
    contracted_edges_mask[idx] = true;

    // 1. CONTRACT EDGE IN THE GRAPH (MIGHT RETURN PARALLEL EDGES DROPPED)
    vector<long> dropped_edges = instance->graph->lemon_contract_edge(idx);

    if (!dropped_edges.empty())
    {
        #ifdef DEBUG
        cout << "edge contraction implied dropping (" << dropped_edges.size()
             << ") parallel edges" << endl;
        #endif

        // fix corresponding model vars at zero
        for (vector<long>::iterator it = dropped_edges.begin();
             it != dropped_edges.end(); ++it)
        {
            fixed_vars.push_back( make_pair(*it, false) );
            model->fix_var(*it, false);
        }
    }

    // 2. FIX VAR = 1 IN THE MODEL (ALSO FIX CONFLICTING EDGES AT 0)
    vector<long> conflicting_vars = model->fix_var(idx, true);

    if (!conflicting_vars.empty())
    {
        #ifdef DEBUG
        cout << "fix var at 1 implied fixing (" << conflicting_vars.size()
             << ") conflicting edges" << endl;
        #endif

        // remove corresponding graph edges
        for (vector<long>::iterator it = conflicting_vars.begin();
             it != conflicting_vars.end(); ++it)
        {
            fixed_vars.push_back( make_pair(*it, false) );
            instance->graph->lemon_delete_edge(*it);
        }
    }
}

IO* LDDA::flush_fixes_to_instance()
{
    /// read current lemon graph and fixed vars to create new IO object

    // TO DO: all

    return 0;
}

void LDDA::print_edge_weights()
{
    // print current weights (in MST subproblem)
    cout << endl << "edge list weights: ";
    for (long idx=0; idx < instance->graph->num_edges; ++idx)
        cout << this->instance->graph->w[idx] << " ";

    cout << endl << "lemonlist weights: ";
    for (long idx=0; idx < instance->graph->num_edges; ++idx)
        cout << (*instance->graph->lemon_weight)[ instance->graph->lemon_edges[idx] ] << " ";
    cout << endl << endl;
}
