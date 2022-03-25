#include <cstdio>
#include <set>
#include <list>
#include <vector>
#include <cmath>
#include <sys/times.h>

#include "ldda_vol.h"

// algorithm setup switches

string VOL_CONFIG_FILE = string("ldda_vol.par");
bool USE_VOL_TIME_LIMIT = true;
double VOL_TIME_LIMIT = 1800;

LDDAVolume::LDDAVolume(IO *instance, KStabModel *model)
: LDDA(instance, model)
{
    this->volume_runtime = -1;
    this->volume_iterations = 0;
    this->volume_bound_log = vector<double>();
    this->volume_multipliers_log = vector< vector<double> >();

    // setup coin-or::vol inner representation of the problem
    this->volp = new VOL_problem(VOL_CONFIG_FILE.c_str());
    volp->psize = 2 * instance->graph->num_edges;  // # primal vars
    volp->dsize = instance->graph->num_edges;      // # dual vars

    // user may override by calling initialize_multipliers()
    this->initial_multipliers_given = false;
}

LDDAVolume::~LDDAVolume()
{
    this->volume_bound_log.clear();
    this->volume_multipliers_log.clear();

    delete volp;
}

bool LDDAVolume::initialize_multipliers(vector<double>& mult)
{
    /// start volume algorithm from the given dual point (lagrangean multipliers)

    cout << "(LDDA-)Vol: Lagrangean Decomposition bound approximation by the Volume Algorithm" << endl << endl;

    if (mult.size() != (unsigned) instance->graph->num_edges)
    {
        cout << endl << "WARNING: WRONG NUMBER OF INITIAL MULTIPLIERS "
             << "(IGNORED AND USED VOLUME DEFAULT)" << endl << endl;

        return false;
    }

    this->initial_multipliers_given = true;

    VOL_dvector& init = volp->dsol;
    init.allocate(volp->dsize);

    cout << "initializing volume at given multipliers" << endl;
    for (int i = 0; i < volp->dsize; ++i)
        init[i] = mult.at(i);

    return true;
}

bool LDDAVolume::run_volume()
{
    if (!this->initial_multipliers_given)
        cout << "(LDDA-)Vol: Lagrangean Decomposition bound approximation by the Volume Algorithm" << endl << endl;

    this->volume_bound = numeric_limits<double>::min();
    this->time_limit_exceeded = false;
    this->start_timer();

    // invoke COIN-OR Vol engine
    long volp_status = volp->solve(*this, this->initial_multipliers_given);
    this->volume_runtime = this->total_time();

    if (volp_status >= 0 || time_limit_exceeded)
    {
        if (time_limit_exceeded)
            cout << "volume time limit exceeded" << endl;
        else
            cout << "volume ended normally" << endl;

        // recompute the violation of the fractional primal solution
        // (adapted from the original UFL example calling COIN-OR Vol)
        vector<double> mismatch = vector<double>(instance->graph->num_edges, 0);
        const VOL_dvector& psol = volp->psol;   // final primal solution

        // mismatch[e] = x[e] - y[e]
        for (long idx=0; idx < instance->graph->num_edges; ++idx)
            mismatch[idx] = psol[idx] - psol[instance->graph->num_edges + idx];

        double violation = 0.0;
        for (long i=0; i < instance->graph->num_edges; ++i)
            violation += fabs(mismatch[i]);
        violation /= (instance->graph->num_edges);
        cout << "average violation of final solution: " << violation << endl;

        // TO DO: if using heuristics, consider additional runs here
        /*
        double heur_val;
        for (long i = 0; i < FINAL_HEURISTIC_RUNS; ++i)
        {
            heur_val = numeric_limits<double>::max();
            this->heuristics(volp, psol, heur_val);
        }
        */

        return true;
    }
    else
    {
        cout << "volume ended prematurely" << endl;
        return false;
    }
}

void inline LDDAVolume::update_edge_weights_if_active(const vector<double>& new_weights)
{
    /***
     * An earlier execution of LDDA might have fixed some elements, so this
     * method is a safe substitute for Graph::update_all_weights(). The 
     * corresponding method in KStabModel is safe as variables simply have
     * their lower and upper bounds set to some fixed value
     */

    for (long idx = 0; idx < this->instance->graph->num_edges; ++idx)
    {
        if ( this->contracted_edges_mask.at(idx) == true || 
             this->removed_edges_mask.at(idx) == true )
        {
            // will save the new weight in vector w but not try to update the
            // lemon weight (which would be raise an error)

            this->instance->graph->w[idx] = new_weights[idx];
        }
        else
            this->instance->graph->update_single_weight(idx, new_weights[idx]);
    }
}

int LDDAVolume::compute_rc(const VOL_dvector& u, VOL_dvector& rc)
{
    /*** 
     * ESSENTIAL USER HOOK (#1 OF 3) IN THE COIN-OR VOL FRAMEWORK
     * Compute "reduced costs": read current multipliers in vector u and
     * determine the objective coefficients of the next lagrangean subproblem
     * (e.g. c_ij - u_j in UFL).
     * See INOC Section 3 for the Lagrangean Decomposition formulation under
     * consideration for finding stable spanning trees.
     */

    // 1. check time limit before starting a new volume iteration

    double current_elapsed_time = total_time();

    if (USE_VOL_TIME_LIMIT && current_elapsed_time > VOL_TIME_LIMIT)
    {
        this->time_limit_exceeded = true;
        return -1;
    }

    cout << " \t\t\t\t\t runtime: " << current_elapsed_time << endl << endl;
    this->volume_iterations++;

    // 2. determine variable coefficients in the objective vector

    const long& m = instance->graph->num_edges;

    for (long idx=0; idx < m; ++idx)
    {
        // mst subproblem (x vars): rc[i] = w_i - u[i] for i in 0...|E|-1
        rc[idx] = this->original_weights[idx] - u[idx];

        // kstab subproblem (y vars): rc[i] = u[i] for i in |E|...2*|E|-1
        rc[m+idx] = u[idx];
    }

    return 0;   // return -1 to abort COIN-OR Vol
}

int LDDAVolume::solve_subproblem( const VOL_dvector& u,    // input: multipliers
                                  const VOL_dvector& rc,   // input: "reduced costs"
                                  double& lcost,           // output: subproblem opt
                                  VOL_dvector& x,          // output: subproblem solution
                                  VOL_dvector& v,          // output: difference between rhs and lhs when plugging x into the relaxed constraints
                                  double& pcost )          // output: objective val of solution in the original problem
{
    /***
     * ESSENTIAL USER HOOK (#2 OF 3) IN THE COIN-OR VOL FRAMEWORK
     * Solve Lagrangean subproblems with "reduced costs" rc.
     * See INOC Section 3 for the Lagrangean Decomposition formulation under
     * consideration for finding stable spanning trees.
     */

    const long& m = instance->graph->num_edges;

    // save current multipliers and prepare new weights for subproblems' objectives
    vector<double> current_multipliers;
    vector<double> mst_weights;
    vector<double> kstab_weights;

    for (long idx=0; idx < m; ++idx)
    {
        current_multipliers.push_back(u[idx]);
        mst_weights.push_back(rc[idx]);
        kstab_weights.push_back(rc[m+idx]);
    }
    volume_multipliers_log.push_back(current_multipliers);

    // 1. SOLVE MST(G, w-lambda^r)

    update_edge_weights_if_active(mst_weights);

    if (instance->graph->mst() == false)
    {
        cout << endl << "UNEXPECTED ERROR: Graph::mst() failed "
             << "- graph not connected?" << endl
             << "infeasible problem instance?" << endl;

        this->problem_solved = true;
        return -1;
    }

    /* For separation of concerns purposes, edges contracted by LDDA (in
     * the LEMON data structure only!) are missing in later spanning trees,
     * so we include them here
     */
    for ( vector<long>::iterator it = contracted_edges.begin();
          it != contracted_edges.end(); ++it )
    {
        instance->graph->mst_weight += mst_weights.at(*it);
        instance->graph->mst_vector[*it] = true;
    }

    // 2. SOLVE KSTAB(\hat{G}, lambda^r)

    model->update_all_weights(kstab_weights);

    if (model->solve(false) <= 0 || model->solution_status != AT_OPTIMUM)
    {
        cout << endl << "UNEXPECTED ERROR: KStabModel::solve() failed" << endl;
        if (model->solution_status == IS_INFEASIBLE)
        {
            cout << "infeasible problem instance" << endl;
            this->problem_solved = true;
        }
        else
            cout << "kstab time limit exceeded?" << endl;

        return -1;
    }

    // 3. SAVE SUBPROBLEM SOLUTIONS IN SINGLE VOL_dvector x

    for (long idx=0; idx < m; ++idx)
    {
        x[idx] = instance->graph->mst_vector.at(idx);
        x[m+idx] = model->solution_vector.at(idx);
    }

    // 4. DETERMINE LAGRANGEAN OBJ (LCOST) AND ORIGINAL OBJ (PCOST)

    lcost = instance->graph->mst_weight + model->solution_weight;
    volume_bound_log.push_back(lcost);
    if (lcost > this->volume_bound)
        this->volume_bound = lcost;

    // TO DO: could not decide whether "the primal solution" should be the x or
    // y solution; preliminary tests gave the exact same behavior
    pcost = 0;
    for (long idx=0; idx < m; ++idx)
        if (model->solution_vector.at(idx) == true)
            pcost += this->original_weights[idx];

    // 5. DETERMINE MISMATCH VECTOR V =  "0 - (X-Y)" FROM THE SUBPROBLEM SOLUTIONS

    for (long idx=0; idx < m; ++idx)
    {
        v[idx] = model->solution_vector.at(idx) - instance->graph->mst_vector.at(idx);
    }

   return 0;   // return -1 to abort COIN-OR Vol
}


// IN:  fractional primal solution (x),
//      best feasible soln value so far (icost)
// OUT: integral primal solution (this->ix) if better than best so far
//      and primal value (this->icost)
int LDDAVolume::heuristics( const VOL_problem& p,
            		 const VOL_dvector& x,
                     double& new_icost)
{
    /***
     * ESSENTIAL USER HOOK (#3 OF 3) IN THE COIN-OR VOL FRAMEWORK
     */

    // TO DO: check if solution is integral, and then feasible
    // TO DO: (MAYBE) if solution is not integral, round some of the
    // mismatching vars and check feasibility

    // TO DO: if running heuristics, give maxweightST primal bound to volume

    return 0;   // return -1 to abort COIN-OR Vol
}
