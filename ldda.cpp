#include "ldda.h"

LDDA::LDDA(IO *instance, Model *model)
{
    this->instance = instance;
    this->model = model;
    this->fixed_vars = vector< pair<long,bool> >();

    // initialize multipliers at zero
    multipliers_log = vector< vector<long> >();
    multipliers_log.push_back( vector<long>(instance->graph->num_edges, 0) );
}

LDDA::LDDA(IO *instance, Model *model, vector<long> initial_multipliers)
{
    this->instance = instance;
    this->model = model;
    this->fixed_vars = vector< pair<long,bool> >();

    // use given values as the initial multipliers
    multipliers_log = vector< vector<long> >();
    multipliers_log.push_back(initial_multipliers);
}

LDDA::~LDDA()
{
    this->multipliers_log.clear();
    this->fixed_vars.clear();
}

long LDDA::edge_deletion_bound()
{
    return 0;
}

long LDDA::edge_contraction_bound()
{
    return 0;
}

long LDDA::vertex_deletion_bound()
{
    return 0;
}

long LDDA::vertex_fix_bound()
{
    return 0;
}

long LDDA::dual_ascent(bool steepest_ascent)
{
    /***
     * Main loop
     * The objective function in the Lagrangean Decomposition formulation is
     * z(\lambda) = min {(w - \lambda) x} + min { \lambda y}, with x taken in
     * the set of spanning trees in G and y taken in the set of cardinality
     * |V|-1 stable sets in the conflict graph \hat{G}.
     */

    bool problem_solved = false;
    bool multipliers_updated = true;

    long iter = 0;

    do
    {
        // 1. SOLVE MST(G, w-lambda^r)
        // get solution x^r and cost Z_x^r

        // 2. SOLVE KSTAB(\hat{G}, lambda^r)
        // get solution y^r and cost Z_y^r

        // 3. COLLECT SET OF INDICES OF VARS WHERE THE SOLUTIONS DON'T MATCH

        // 4. TREAT EXCEPTIONS

        // 4.1 THE SET IS EMPTY: THE SOLUTIONS MATCH AND THE PROBLEM IS SOLVED

        // 4.2 IF THE MST SOLUTION IS STABLE, IT IS AN INTEGER FEASIBLE POINT

        // if it has the same lambda^r cost as the kstab, than it is optimal

        // 4.3 IF THE KSTAB SOLUTION IS ACYCLIC, IT IS AN INTEGER FEASIBLE POINT

        // if it has the same lambda^r cost as the kstab, than it is optimal

        // 5. UNLESS WE SOLVED THE PROBLEM, CHOOSE ELEMENT FOR DUAL ASCENT
        // switch between steepest ascent and first ascent
        // switch between x^r_e = 1 and y^r_e = 0 and vice versa
        // determine delta^r_e and del^r_e
        // if probing finds an infeasible problem, fix corresponding variable
        // if the maximum step size along e is zero, pick a different element

        // 6. UPDATE MULTIPLIER \lambda^{r+1}_e, COPY THE OTHERS

        // 7. SCREEN LOG
    }
    while (!problem_solved && multipliers_updated)

    if (problem_solved)
    {
        // original problem is solved
    }
    else
    {
        // original problem not solved, but no ascent was possible
    }

    //instance->graph->mst();   // ignoring return value, assuming a spanning tree exists
    //model->solve(true);

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
