#include "mstcc_model.h"
#include "mstcc_cut_generator.h"

#define ORTHOGONALITY_TOL_IN_LPR 0.1
#define SEC_VIOLATION_TOL_IN_LPR 0.0001
#define SEC_SEPARATION_PRECISION_IN_LPR 14   // <= 14 without changing everything to long double (which gurobi ignores)

// algorithm setup switches (set the first 3 options to true to compute the LP relaxation faster)
bool USE_DEGREE_CONSTRAINTS_A_PRIORI = false;
bool USE_FAST_FOLKLORE_CUT_IN_LPR = false;
bool USE_FAST_INTEGER_CUT_IN_LPR = true;

double TIME_LIMIT_IN_IP = 3600;

StableSpanningTreeModel::StableSpanningTreeModel(IO *instance)
: KStabModel(instance)
{
    this->lp_bound = this->lp_runtime = -1;

    if (USE_DEGREE_CONSTRAINTS_A_PRIORI)
    {
        // enforce a priori that all vertices should have degree >= 1
        // (implied by SEC), to make LP relaxation faster
        ostringstream cname;
        for (long u=0; u < instance->graph->num_vertices; u++)
        {
            GRBLinExpr cut_edges = 0;
            for (list<long>::iterator it = instance->graph->adj_list[u].begin();
                it != instance->graph->adj_list[u].end(); ++it)
            {
                long edge_idx = instance->graph->index_matrix[u][*it];
                cut_edges += x[edge_idx];
            }

            cname.str("");
            cname << "MSTCC_degree_" << u;
            model->addConstr(cut_edges >= 1, cname.str());
        }
        model->update();
    }
}

StableSpanningTreeModel::~StableSpanningTreeModel()
{
    /// NOTHING HERE
}

int StableSpanningTreeModel::solve(bool logging)
{
    /***
     * Polymorphic override of the kstab solve() method.
     * Finds a min weight solution of the natural IP formulation for MSTCC:
     * kstab in the conflict graph + subtour elimination constraints (SEC) in
     * the original one. Returns the number of solutions found, while further
     * information is stored in object fields.
     */

    try
    {
        /*
        // turn off all built-in cut generators?
        model->set(GRB_IntParam_Cuts, 0);

        // turn off all preprocessing and heuristics?
        model->set(GRB_IntParam_Presolve, 0);
        model->set(GRB_IntParam_PrePasses, 0);
        model->set(GRB_DoubleParam_PreSOS1BigM, 0);
        model->set(GRB_DoubleParam_PreSOS2BigM, 0);
        model->set(GRB_IntParam_PreSparsify, 0);
        model->set(GRB_IntParam_PreCrush, 1);
        model->set(GRB_IntParam_DualReductions, 0);
        model->set(GRB_IntParam_Aggregate, 0);

        model->set(GRB_DoubleParam_Heuristics, 0);
        */

        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

        model->set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_IN_IP);

        // should disable presolve reductions that affect user cuts
        model->set(GRB_IntParam_PreCrush, 1);

        // must set parameter indicating presence of lazy constraints
        model->set(GRB_IntParam_LazyConstraints, 1);

        // set callback to separate SECs and solve IP
        SSTCutGenerator cutgen = SSTCutGenerator(model, x, instance);
        model->setCallback(&cutgen);
        model->optimize();

        return this->save_optimization_status();
    }
    catch(GRBException e)
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return 0;
    }
    catch(...)
    {
        cout << "Unexpected error during optimization inside StableSpanningTreeModel::solve()" << endl;
        return 0;
    }

}

bool StableSpanningTreeModel::solve_lp_relax(bool logging)
{
    /***
     * Solves the lp relaxation of the natural IP formulation for MSTCC:
     * kstab in the conflict graph + subtour elimination constraints (SEC) in
     * the original one. Returns true if the bound was computed, and false if
     * the LP formulation is already infeasible.
     */

    try
    {
        // turn off all gurobi cut generators
        model->set(GRB_IntParam_Cuts, 0);

        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

        // this cut generator object is used to find violated SECs only (not via gurobi's callback)
        SSTCutGenerator cutgen = SSTCutGenerator(model, x, instance);

        // make vars continuous
        for (long i=0; i < instance->graph->num_edges; i++)
            x[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);

        // calculating wall clock time to solve LP relaxation
        struct timeval *clock_start = (struct timeval *) malloc(sizeof(struct timeval));
        struct timeval *clock_stop  = (struct timeval *) malloc(sizeof(struct timeval));
        gettimeofday(clock_start, 0);

        // solve LP to optimality; then iterate reoptimizing and separating SECs
        model->optimize();
        bool model_updated = true;
        this->lp_passes = 1;

        while (model_updated && model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            #ifdef DEBUG_LPR
            cout << "LP relaxation pass #" << lp_passes << " (bound = "
                 << model->get(GRB_DoubleAttr_ObjVal) << ")" << endl;
            #endif

            model_updated = cutgen.separate_lpr();

            if (model_updated)
            {
                // reoptimize
                model->update();
                model->optimize();
                this->lp_passes++;
            }
        }

        // LP relaxation time
        gettimeofday(clock_stop, 0);
        unsigned long clock_time = 1.e6 * (clock_stop->tv_sec - clock_start->tv_sec) +
                                          (clock_stop->tv_usec - clock_start->tv_usec);
        this->lp_runtime = ((double)clock_time / (double)1.e6);
        free(clock_start);
        free(clock_stop);

        // loop might have broken because no violated SEC exists or because the problem became infeasible
        if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            this->lp_bound = model->get(GRB_DoubleAttr_ObjVal);
            
            // restore IP model
            model->set(GRB_IntParam_Cuts, -1);
            for (long i=0; i < instance->graph->num_edges; i++)
                x[i].set(GRB_CharAttr_VType, GRB_BINARY);

            return true;
        }
        else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
        {
            cout << "LP relaxation infeasible!" << endl;
            cout << "Model runtime: " << lp_runtime << endl;
            return false;
        }
        else
        {
            cout << "Unexpected error: solve_lp_relax() got neither optimal nor infeasible model" << endl;
            cout << "Model runtime: " << lp_runtime << endl;
            return false;
        }
    }
    catch(GRBException e)
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return false;
    }
}
