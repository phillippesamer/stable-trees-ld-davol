#include "kstab_model.h"

KStabModel::KStabModel(IO *instance)
{
    this->instance = instance;
    this->fixed_cardinality = instance->graph->num_vertices-1;

    try
    {
        this->env = new GRBEnv();
        this->model = new GRBModel(*env);

        create_variables();
        create_constraints();
        create_objective();

        //model->write("kstab.lp");
    }
    catch(GRBException e)
    {
        cout << "Model construction error, code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}

KStabModel::~KStabModel()
{
    delete[] x;
    delete model;
    delete env;
}

void KStabModel::create_variables()
{
    char buffer[50];

    // binary selection of each vertex (edge in the original graph)
    x = new GRBVar[instance->graph->num_edges];
    for (long i=0; i < instance->graph->num_edges; i++)
    {
        sprintf(buffer, "x_%ld", i);
        x[i] = model->addVar(0.0, 1.0, 1.0, GRB_BINARY, buffer);
    }
    model->update();
}

void KStabModel::create_constraints()
{
    ostringstream cname;

    // 1. FIXED CARDINALITY CONSTRAINT
    GRBLinExpr fix_card = 0;
    for (long i=0; i < instance->graph->num_edges; i++)
        fix_card += x[i];

    cname.str("");
    cname << "C1_fixed_cardinality";
    model->addConstr(fix_card == this->fixed_cardinality, cname.str());

    // 2. EDGE INEQUALITIES
    for (long e=0; e < instance->num_conflicts; e++)
    {
        long u = instance->conflicts[e].first;
        long v = instance->conflicts[e].second;

        GRBLinExpr edge_ineq = 0;
        edge_ineq += x[u];
        edge_ineq += x[v];

        cname.str("");
        cname << "C2_edge_" << u << "_" << v;
        model->addConstr(edge_ineq <= 1, cname.str());
    }

    model->update();
}

void KStabModel::create_objective()
{
    GRBLinExpr objective_expression = 0;

    for (long i = 0; i < instance->graph->num_edges; i++)
        objective_expression += (instance->graph->w[i])*x[i];

    model->setObjective(objective_expression, GRB_MINIMIZE);
    model->update();
}

int KStabModel::solve(bool logging)
{
    try
    {
        // turn off all gurobi cut generators
        //model->set(GRB_IntParam_Cuts, 0);

        // turn off gurobi presolve and heuristics
        //model->set(GRB_IntParam_Presolve, 0);
        //model->set(GRB_IntParam_PrePasses, 0);
        //model->set(GRB_DoubleParam_Heuristics, 0);

        // need to disable presolve reductions that affect user cuts
        //model->set(GRB_IntParam_PreCrush, 1);

        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

        model->optimize();
        this->solution_runtime = model->get(GRB_DoubleAttr_Runtime);

        // store results in this object
        if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        {
            this->solution_status = AT_OPTIMUM;
            this->solution_vector.clear();

            // the weights are integral, so the optimal value must also be
            double optval = model->get(GRB_DoubleAttr_ObjVal);
            this->solution_weight = (long) optval;
            this->solution_dualbound = (long) optval;
            if (floor(optval) - optval != 0)
            {
                this->solution_weight = std::floor(optval + ZERO_TOL);
                this->solution_dualbound = std::floor(optval + ZERO_TOL);
                cout << "NB: model optimal value is " << optval
                     << " (not integral!); rounded to " << solution_weight
                     << endl;
            }

            // save bool vector of this solution

            /*
            #ifdef DEBUG
                cout << "Model optimal vector: " << endl;
            #endif
            */

            for (long i=0; i < instance->graph->num_edges; ++i)
            {
                // NB: gurobi vars are floating point numbers, allowing +0 and -0!
                // The vector<bool> below is safer and easier to use

                if (this->x[i].get(GRB_DoubleAttr_X) > ZERO_TOL)
                    this->solution_vector.push_back(true);
                else
                    this->solution_vector.push_back(false);

                /*
                #ifdef DEBUG
                    cout << this->solution_vector.back();
                #endif
                */
            }

            // returning number of solutions (likely sub-optimal), in case it's useful later
            return model->get(GRB_IntAttr_SolCount);
        }
        else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
        {
            this->solution_status = IS_INFEASIBLE;

            this->solution_weight = numeric_limits<long>::max();
            this->solution_dualbound = numeric_limits<long>::max();
            this->solution_vector.clear();

            cout << "Model infeasible! (runtime "
                 << solution_runtime << ")" << endl;

            return 0;
        }
        else if (model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
        {
            this->solution_status = STATUS_UNKNOWN;

            this->solution_weight = model->get(GRB_DoubleAttr_ObjVal);
            this->solution_dualbound = model->get(GRB_DoubleAttr_ObjBound);

            cout << "Time limit exceeded (" << solution_runtime << ")" << endl;
            cout << "Dual bound " << model->get(GRB_DoubleAttr_ObjBound) 
                 << ", primal bound " << model->get(GRB_DoubleAttr_ObjVal) 
                 << " (MIP gap " << model->get(GRB_DoubleAttr_MIPGap) << "%)" 
                 << endl;

            return 0;
        }
        else
        {
            this->solution_status = STATUS_UNKNOWN;
            this->solution_vector.clear();

            cout << "Unexpected error: solve() got neither optimal nor infeasible model" << endl;

            return 0;
        }
    }
    catch(GRBException e)
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return 0;
    }
}

double KStabModel::runtime()
{
    return model->get(GRB_DoubleAttr_Runtime);
}


void KStabModel::set_time_limit(double tl)
{
    model->set(GRB_DoubleParam_TimeLimit, tl);
}

void KStabModel::update_single_weight(long idx, long new_weight)
{
    x[idx].set(GRB_DoubleAttr_Obj, new_weight);
    model->update();

    /*
    #ifdef DEBUG
        model->write("kstab_2.lp");
    #endif
    */
}

void KStabModel::update_all_weights(vector<long> new_weights)
{
    GRBLinExpr objective_expression = 0;

    for(long i = 0; i < instance->graph->num_edges; i++)
        objective_expression += (new_weights[i])*x[i];

    model->setObjective(objective_expression, GRB_MINIMIZE);
    model->update();

    /*
    #ifdef DEBUG
        model->write("kstab_3.lp");
    #endif
    */
}

pair<ModelStatus,long> KStabModel::probe_var(long probe_idx, bool probe_value)
{
    // change given var to continuous and update lb=ub=probe
    this->fix_var(probe_idx, probe_value);

    /*
    #ifdef DEBUG
        model->write("kstab_2.lp");
    #endif
    */

    // TO DO: set output 0 and then 1 again? 
    //model->set(GRB_IntParam_OutputFlag, 0);
    model->optimize();
    //model->set(GRB_IntParam_OutputFlag, 1);

    // save optimization status (optimal/infeasible) and result, if optimal
    ModelStatus status = STATUS_UNKNOWN;
    long result = numeric_limits<long>::max();

    if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        status = AT_OPTIMUM;

        // the weights are integral, so the optimal value must also be
        double optval = model->get(GRB_DoubleAttr_ObjVal);
        result = (long) optval;
        if (floor(optval) - optval != 0)
        {
            result = std::floor(optval + ZERO_TOL);

            cout << "NB: probe_var(idx=" << probe_idx << ", value=" << probe_value << ") ";
            cout << "got optimal value " << optval << " (not integral!); rounded to "
                 << result << endl;
        }
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        status = IS_INFEASIBLE;
        result = numeric_limits<long>::max();
    }
    else
        cout << "Unexpected error: probe_var(idx=" << probe_idx << ", value=" << probe_value << ") got neither optimal nor infeasible model" << endl;

    // restore all variables: binary, lb=0, ub=1
    x[probe_idx].set(GRB_DoubleAttr_LB, 0);
    x[probe_idx].set(GRB_DoubleAttr_UB, 1);
    x[probe_idx].set(GRB_CharAttr_VType, GRB_BINARY);

    if (probe_value == true)
    {
        for (list<long>::iterator it = instance->conflict_graph_adj_list[probe_idx].begin(); 
             it != instance->conflict_graph_adj_list[probe_idx].end(); ++it)
        {
            x[*it].set(GRB_DoubleAttr_LB, 0);
            x[*it].set(GRB_DoubleAttr_UB, 1);
            x[*it].set(GRB_CharAttr_VType, GRB_BINARY);
        }
    }

    model->update();

    /*
    #ifdef DEBUG
        model->write("kstab_3.lp");
    #endif
    */

    return make_pair(status,result);
}

vector<long> KStabModel::fix_var(long fix_idx, bool fix_value)
{
    /// change given var to continuous and update lb=ub=fix_value
    double probe_dbl = fix_value == true ? 1.0 : 0.0;         // redundant, but...
    x[fix_idx].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    x[fix_idx].set(GRB_DoubleAttr_LB, probe_dbl);
    x[fix_idx].set(GRB_DoubleAttr_UB, probe_dbl);

    /*
    #ifdef DEBUG
        cout << "fixed x[" << fix_idx << "] = " << fix_value << endl;
    #endif
    */

    vector<long> conflicting_vars = vector<long>();

    // if fixing at 1, change vars corresponding to neighbours: continuous, lb=ub=0
    if (fix_value == true)
    {
        for (list<long>::iterator it = instance->conflict_graph_adj_list[fix_idx].begin(); 
             it != instance->conflict_graph_adj_list[fix_idx].end(); ++it)
        {
            x[*it].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
            x[*it].set(GRB_DoubleAttr_LB, 0);
            x[*it].set(GRB_DoubleAttr_UB, 0);

            /*
            #ifdef DEBUG
                cout << "fixed x[" << fix_idx << "] at zero" << endl;
            #endif
            */

            conflicting_vars.push_back(*it);
        }
    }

    model->update();

    return conflicting_vars;
}
