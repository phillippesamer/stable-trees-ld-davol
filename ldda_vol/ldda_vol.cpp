#include <cstdio>
#include <set>
#include <list>
#include <vector>
#include <cmath>
#include <sys/times.h>

#include "ldda_vol.h"

LDDAVolume::LDDAVolume(IO *instance, KStabModel *model)
: LDDA(instance, model)
{
    // might be overridden by config file
    this->initialized_mult = false;
    this->h_iter = 10;
    this->volume_precision = 14;

    this->volume_runtime = -1;
    this->volume_iterations = 0;
    this->volume_bound_log = vector<double>();
    this->volume_multipliers_log = vector< vector<double> >();

    this->volp = NULL;
}

LDDAVolume::~LDDAVolume()
{
    this->volume_bound_log.clear();
    this->volume_multipliers_log.clear();

    if (volp != NULL)
        delete volp;
}

bool LDDAVolume::read_volume_config_file(string filename)
{
    // original IO method in the UFL example file in COIN-OR Vol

    // TO DO: update code for readability to use c++ api

    this->config_file = filename;

    // 1. OPTIONAL USER PARAMETERS

    char s[500];
    FILE * file = fopen(filename.c_str(), "r");
    if (!file)
    {
        cout << "could not open volume config file: " << filename << endl;
        return false;
    }

    while (fgets(s, 500, file))
    {
        unsigned long len = strlen(s) - 1;
        if (s[len] == '\n')
            s[len] = 0;
        string ss = s;

        if (ss.find("dualfile") == 0)
        {
            unsigned long j = ss.find("=");
            long j1 = ss.length() - j + 1;
            dualfile = ss.substr(j+1, j1);
        }
        else if (ss.find("dual_savefile") == 0)
        {
            unsigned long j = ss.find("=");
            long j1 = ss.length() - j + 1;
            dual_savefile = ss.substr(j+1, j1);
        }
        else if (ss.find("h_iter") == 0)
        {
            unsigned long i = ss.find("=");  
            h_iter = atol(&s[i+1]);
        }
        else if (ss.find("volume_precision") == 0)
        {
            unsigned long i = ss.find("=");  
            volume_precision = atoi(&s[i+1]);
        }
    }

    fclose(file);

    // 2. PREPARE COIN-OR VOL INNER REPRESENTATION OF THE PROBLEM

    this->volp = new VOL_problem(filename.c_str());
    volp->psize = 2 * this->instance->graph->num_edges;  // # primal vars
    volp->dsize = this->instance->graph->num_edges;      // # dual vars

    // check if the user defined initial multipliers
    // TO DO: replace this by a new constructor?
    if (this->dualfile.length() > 0)
    {
        initialized_mult = true;

        // read dual solution
        VOL_dvector& dinit = volp->dsol;
        dinit.allocate(volp->dsize);

        FILE *file = fopen(this->dualfile.c_str(), "r");
        if (!file)
        {
            cout << "could not open initial multipliers file: " << this->dualfile.c_str() << endl;
            return false;
        }

        int idummy;
        cout << "initializing volume at:" << endl;
        for (int i = 0; i < volp->dsize; ++i)
        {
            fscanf(file, "%i%lf", &idummy, &dinit[i]);
            cout << "lambda[" << i << "] = " << dinit[i] << ", ";
        }
        cout << endl;

        fclose(file);
    }

    return true;
}

bool LDDAVolume::run_volume()
{
    // adapted from the original UFL example calling COIN-OR Vol

    // TO DO: give maxweightST primal bound to volume
    // TO DO: set correct filename of the initial multipliers on .par file

    this->start_timer();

    // invoke COIN-OR Vol engine
    if (volp->solve(*this, this->initialized_mult) < 0)
    {
        cout << "COIN-OR Vol::solve() failed" << endl;
        this->volume_runtime = total_time();
        return false;
    }
    else
    {
        // recompute the violation of the fractional primal solution
        vector<double> mismatch = vector<double>(instance->graph->num_edges, 0);
        const VOL_dvector& psol = volp->psol;   // final primal solution

        for (long idx=0; idx < instance->graph->num_edges; ++idx)
        {
            // mismatch[e] = x[e] - y[e]
            mismatch[idx] = psol[idx] - psol[instance->graph->num_edges + idx];
        }

        double violation = 0.0;
        for (long i=0; i < instance->graph->num_edges; ++i)
            violation += fabs(mismatch[i]);
        violation /= (instance->graph->num_edges);
        cout << "average violation of final solution: " << violation << endl;

        if (this->dual_savefile.length() > 0)
        {
            // save dual solution
            FILE* file = fopen(this->dual_savefile.c_str(), "w");
            const VOL_dvector& u = volp->dsol;
            for (int i = 0; i < u.size(); ++i)
                fprintf(file, "%8i  %f\n", i+1, u[i]);
            fclose(file);
        }

        // TO DO: additional runs of the heuristics?
        /*
        double heur_val;
        for (long i = 0; i < this->h_iter; ++i)
        {
            heur_val = DBL_MAX;
            this->heuristics(volp, psol, heur_val);
        }
        */
    }

    this->volume_runtime = this->total_time();
    cout << "volume bound: " << volp->value
         << " (runtime " << fixed << volume_runtime << ")" << endl;

    return true;
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
                                  double& pcost )          // output: objective val of x in the original problem
{
    /***
     * ESSENTIAL USER HOOK (#2 OF 3) IN THE COIN-OR VOL FRAMEWORK
     * Solve Lagrangean subproblems with "reduced costs" rc.
     * See INOC Section 3 for the Lagrangean Decomposition formulation under
     * consideration for finding stable spanning trees.
     */

    this->volume_iterations++;

    const long& m = instance->graph->num_edges;

    // save current multipliers and new weights for the objectives of the subproblems
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

    // 1. SOLVE MST SUBPROBLEM WITH WEIGHTS rc

    //instance->graph->update_all_weights(mst_weights);

    // need to account for contracted/deleted edges


    // 2. SOLVE KSTAB SUBPROBLEM WITH WEIGHTS rc

    //model->update_all_weights(kstab_weights);

    // need to account for contracted/deleted edges



    // 3. SAVE SUBPROBLEM SOLUTIONS IN SINGLE VOL_dvector x  
    // need to account for contracted/deleted edges



    // 4. DETERMINE LAGRANGEAN OBJ (LCOST) AND ORIGINAL OBJ (PCOST)
    // need to account for contracted/deleted edges
    // volume_bound_log.push_back(lcost?);



    // 5. DETERMINE MISMATCH VECTOR V =  "0 -(X-Y)" WITH THE SUBPROBLEM SOLUTIONS
    // need to account for contracted/deleted edges



    /*
   int i,j;

   lcost = 0.0;
   for (i = 0; i < ncust; ++i) {
      lcost += u[i];
      v[i]=1;
   }

   VOL_ivector sol(nloc + nloc*ncust);

   // produce a primal solution of the relaxed problem
   const double * rdist = rc.v + nloc;
   double sum;
   int k=0, k1=0;
   double value=0.;
   int xi;
   for ( i=0; i < nloc; ++i ) {
     sum=0.;
     for ( j=0; j < ncust; ++j ) {
       if ( rdist[k]<0. ) sum+=rdist[k];
       ++k;
     }
     if (fix[i]==0) xi=0;
     else 
       if (fix[i]==1) xi=1;
       else 
	 if ( fcost[i]+sum >= 0. ) xi=0;
	 else xi=1;
     sol[i]=xi;
     value+=(fcost[i]+sum)*xi;
     for ( j=0; j < ncust; ++j ) {
       if ( rdist[k1] < 0. ) sol[nloc+k1]=xi;
       else sol[nloc+k1]=0;
       ++k1;
     }
   }

   lcost += value;

   pcost = 0.0;
   x = 0.0;
   for (i = 0; i < nloc; ++i) {
     pcost += fcost[i] * sol[i];
     x[i] = sol[i];
   }

   k = 0;
   for ( i=0; i < nloc; i++){
     for ( j=0; j < ncust; j++){
       x[nloc+k]=sol[nloc+k];
       pcost+= dist[k]*sol[nloc+k];
       v[j]-=sol[nloc+k];
       ++k;
     }
   }
   */

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

    return 0;   // return -1 to abort COIN-OR Vol
}
