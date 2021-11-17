#include "ldda.h"

LDDA::LDDA(IO *instance, Model *model)
{
    this->instance = instance;
    this->model = model;
    this->fixed_vars = vector< pair<long,bool> >();
    this->bound_log = vector<long>();
    this->solution_pool = vector< pair<long, vector<bool> > >();

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

    // initialize weights with lagrangean multipliers

    // args: source1 begin, source1 end, source2 begin, dest begin, operator
    transform( instance->graph->w.begin(),
               instance->graph->w.end(),
               multipliers_log[0].begin(),
               instance->graph->w.begin(),
               minus<long>() );

    #ifdef DEBUG
        cout << "w contains: ";
        for (vector<long>::iterator it = instance->graph->w.begin(); it != instance->graph->w.end(); ++it)
            cout << *it << " ";
        cout << endl;
    #endif

    instance->graph->update_all_weights(instance->graph->w);

    model->update_all_weights(multipliers_log[0]);

    // TO DO: ANYTHING ELSE ON THE FIRST RUN?

    bool problem_solved = false;
    bool multipliers_updated = true;
    long iter = 0;

    do
    {
        // 1. SOLVE MST(G, w-lambda^r)
        if (instance->graph->mst() == false)
        {
            cout << "Graph::mst() failed (graph is not connected)" << endl;
            cout << "infeasible problem instance" << endl;
            return false;
        }
        //long mst_weight = instance->graph->mst_weight;
        //vector<bool> mst_vector = instance->graph->mst_vector;
        //double mst_runtime = instance->graph->mst_runtime;

        // 2. SOLVE KSTAB(\hat{G}, lambda^r)
        if (model->solve(true) <= 0 || model->solution_status != AT_OPTIMUM)
        {
            cout << "Model::solve() failed" << endl;
            cout << "infeasible problem instance" << endl;
            return false;
        }
        //long kstab_weight = model->solution_weight;
        //vector<bool> kstab_vector = model->solution_vector;
        //double kstab_runtime = model->solution_runtime;

        // 3. COLLECT SET OF INDICES OF VARS WHERE THE SOLUTIONS DON'T MATCH
        vector<long> mismatch;
        for (long idx=0; idx < instance->graph->num_edges; ++idx)
        {
            if (instance->graph->mst_vector[idx] != model->solution_vector[idx])
                mismatch.push_back(idx);
        }

        // 4. TREAT EXCEPTIONS

        // 4.1 THE SET IS EMPTY: THE SOLUTIONS MATCH AND THE PROBLEM IS SOLVED
        if (mismatch.empty())
        {
            cout << "the solutions to both subproblems are the same" << endl;
                 << "problem solved to optimality" << endl;

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
                     << "problem solved to optimality" << endl;

                return true;
            }
        }

        // 4.3 IF THE KSTAB SOLUTION IS ACYCLIC, IT IS AN INTEGER FEASIBLE POINT
        if ( instance->test_acyclic(model->solution_vector) )
        {
            // store integer feasible solution
            
            // if it has the same (w-lambda^r) cost as the mst, than it is optimal
        }

        // 5. UNLESS WE SOLVED THE PROBLEM, CHOOSE ELEMENT FOR DUAL ASCENT
        // switch between steepest ascent and first ascent
        // switch between x^r_e = 1 and y^r_e = 0 and vice versa
        // determine delta^r_e and del^r_e
        // if probing finds an infeasible problem, fix corresponding variable (save this on fixed_vars), pick new one
        // if the maximum step size along e is zero, pick a different element

        // 6. UPDATE MULTIPLIER \lambda^{r+1}_e, COPY THE OTHERS

        // 7. SCREEN LOG
    }
    while (!problem_solved && multipliers_updated);

    if (problem_solved)
    {
        // original problem is solved
    }
    else
    {
        // original problem not solved, but no ascent was possible
    }





    // mst probing
    /*
    cout << endl;
    for (long i=0; i < instance->num_edges(); i++)
    {
        pair<bool,long> probing = instance->graph->mst_probing_var(i,false);
        cout << "probing x[" << i << "]=" << 0 << " gives (" << probing.first << "," << probing.second << ")" << endl;
    }
    cout << endl;
    */

    cout << endl;
    for (long i=0; i < instance->num_edges(); i++)
    {
        pair<ModelStatus,double> probing = model->probe_var(i,true);
        if (probing.first == AT_OPTIMUM)
        {
            cout << "probing var x[" << i << "] = " << true << " gives ObjVal=" << probing.second 
            << " (runtime: " << model->runtime() << " s)" << endl;
        }
        else if (probing.first == IS_INFEASIBLE)
        {
            cout << "probing var x[" << i << "] = 0 gives an infeasible model" 
            << " (runtime: " << model->runtime() << " s)" << endl;
        }
    }
    cout << endl;

    /*
    long m = instance->graph->num_edges;

    print_edge_weights();

    for (long i=0; i<m; ++i)
        instance->graph->update_single_weight(i,10*i);
    print_edge_weights();

    instance->graph->update_all_weights( vector<long>(m, -1) );
    print_edge_weights();

    */

    return true;
}

long LDDA::edge_deletion_bound()
{
    // TO DO: all
    
    return 0;
}

long LDDA::edge_contraction_bound()
{
    // TO DO: all
    
    return 0;
}

long LDDA::vertex_deletion_bound()
{
    // TO DO: all
    
    return 0;
}

long LDDA::vertex_fix_bound()
{
    // TO DO: all

    return 0;
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
