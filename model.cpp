#include "model.h"

Model::Model(IO *instance)
{
    this->instance = instance;
    this->fixed_cardinality = instance->graph->num_vertices-1;

    try
    {
        this->env = new GRBEnv();
        this->model = new GRBModel(*env);

        cout << "creating variables..." << endl;
        create_variables();

        cout << "creating constraints..." << endl;
        create_constraints();

        cout << "setting objective function..." << endl;
        create_objective();

        //model->write("kstab.lp");
    }
    catch(GRBException e)
    {
        cout << "Model construction error, code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
}

Model::~Model()
{
    delete[] x;
    delete model;
    delete env;
}

void Model::create_variables()
{
    char buffer[50];

    // binary selection of each vertex (edge in the original graph)
    x = new GRBVar[instance->graph->num_edges];
    for(long i=0; i < instance->graph->num_edges; i++)
    {
        sprintf(buffer, "x_%ld", i);
        x[i] = model->addVar(0.0, 1.0, 1.0, GRB_BINARY, buffer);
    }
    model->update();
}

void Model::create_constraints()
{
    ostringstream cname;

    // 1. FIXED CARDINALITY CONSTRAINT
    GRBLinExpr fix_card = 0;
    for(long i=0; i < instance->graph->num_edges; i++)
        fix_card += x[i];

    cname.str("");
    cname << "C1_fixed_cardinality";
    model->addConstr(fix_card == this->fixed_cardinality, cname.str());

    // 2. EDGE INEQUALITIES
    for(long e=0; e < instance->num_conflicts; e++)
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

void Model::create_objective()
{
    GRBLinExpr objective_expression = 0;

    for(long i = 0; i < instance->graph->num_edges; i++)
        objective_expression += (instance->graph->w[i])*x[i];

    model->setObjective(objective_expression, GRB_MINIMIZE);
    model->update();
}

int Model::solve()
{
    try
    {
        // turn off all gurobi cut generators
        //model->getEnv().set(GRB_IntParam_Cuts, 0);

        // turn off gurobi presolve and heuristics
        //model->getEnv().set(GRB_IntParam_Presolve, 0);
        //model->getEnv().set(GRB_IntParam_PrePasses, 0);
        //model->getEnv().set(GRB_DoubleParam_Heuristics, 0);

        // need to disable presolve reductions that affect user cuts
        //model->getEnv().set(GRB_IntParam_PreCrush, 1);

        model->optimize();
        return model->get(GRB_IntAttr_SolCount);
    }
    catch(GRBException e)
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        return 0;
    }
}

double Model::runtime()
{
    return model->get(GRB_DoubleAttr_Runtime);
}

pair<ProbeStatus,double> Model::probe_var_at_zero(long var_idx)
{
    // change var: continuous, lb=ub=0
    x[var_idx].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    x[var_idx].set(GRB_DoubleAttr_LB, 0);
    x[var_idx].set(GRB_DoubleAttr_UB, 0);

    model->update();

    /*
    #ifdef DEBUG
        cout << endl << endl << "x[" << var_idx << "] now set to VType=" << 
            x[var_idx].get(GRB_CharAttr_VType) << " (see .lp file!)" << endl;
        model->write("kstab_2.lp");
    #endif
    */

    model->set(GRB_IntParam_OutputFlag, 0);
    model->optimize();
    model->set(GRB_IntParam_OutputFlag, 1);

    // save optimization status (optimal/infeasible) and result, if optimal
    ProbeStatus status = PROBE_UNKNOWN;
    double result = numeric_limits<double>::max();

    if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        status = PROBE_OPTIMAL;
        result = model->get(GRB_DoubleAttr_ObjVal);
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        status = PROBE_INFEASIBLE;
        result = numeric_limits<double>::max();
    }
    else
        cout << "Unexpected error: probe_var_at_zero got neither optimal nor infeasible model" << endl;

    // restore original var: binary, lb=0, ub=1
    x[var_idx].set(GRB_DoubleAttr_LB, 0);
    x[var_idx].set(GRB_DoubleAttr_UB, 1);
    x[var_idx].set(GRB_CharAttr_VType, GRB_BINARY);

    model->update();

    /*
    #ifdef DEBUG
        cout << endl << "x[" << var_idx << "] set back to VType=" << 
            x[var_idx].get(GRB_CharAttr_VType) << " (see .lp file!)" << endl;
        model->write("kstab_3.lp");
    #endif
    */

    return make_pair(status,result);
}

pair<ProbeStatus,double> Model::probe_var_at_one(long var_idx)
{
    // change given var: continuous, lb=ub=1
    x[var_idx].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    x[var_idx].set(GRB_DoubleAttr_LB, 1);
    x[var_idx].set(GRB_DoubleAttr_UB, 1);

    // change vars corresponding to adjacent vertices: continuous, lb=ub=0
    cout << "just looking at the conflict graph adjacency list: ";
    cout << instance->conflict_graph_adj_list.size() << " entries:" << endl;
    for (long i=0; i < instance->graph->num_edges; ++i)
    {
        cout << "list[" << i << "] = { ";
        for (list<long>::iterator it = instance->conflict_graph_adj_list[i].begin(); 
             it != instance->conflict_graph_adj_list[i].end(); ++it)
        {
            cout << *it << " ";
        }
        cout << "}" << endl;
    }

    model->update();

    /*
    #ifdef DEBUG
        cout << endl << endl << "x[" << var_idx << "] now set to VType=" << 
            x[var_idx].get(GRB_CharAttr_VType) << " (see .lp file!)" << endl;
        model->write("kstab_2.lp");
    #endif
    */

    //model->set(GRB_IntParam_OutputFlag, 0);
    //model->optimize();
    //model->set(GRB_IntParam_OutputFlag, 1);
//////////////////////

    // save optimization status (optimal/infeasible) and result, if optimal
    ProbeStatus status = PROBE_UNKNOWN;
    double result = numeric_limits<double>::max();

    if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        status = PROBE_OPTIMAL;
        result = model->get(GRB_DoubleAttr_ObjVal);
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        status = PROBE_INFEASIBLE;
        result = numeric_limits<double>::max();
    }
    else
        cout << "Unexpected error: probe_var_at_zero got neither optimal nor infeasible model" << endl;

    // restore all variables: binary, lb=0, ub=1
    // 

    return make_pair(status,result);
}
