#include "model.h"

void Model::create_variables()
{
    char buffer[50];

    // binary selection of each vertex
    x = new GRBVar[instance->num_vertices];
    for(int i=0; i < instance->num_vertices; i++)
    {
        sprintf(buffer, "x_%d", i);
        x[i] = model->addVar(0.0, 1.0, 1.0, GRB_BINARY, buffer);
    }
    model->update();
}

void Model::create_constraints()
{
    ostringstream cname;

    // 1. FIXED CARDINALITY CONSTRAINT
    GRBLinExpr fixed_cardinality = 0;
    for(int i=0; i < instance->num_vertices; i++)
        fixed_cardinality += x[i];

    cname.str("");
    cname << "C1_fixed_cardinality";
    model->addConstr(fixed_cardinality == card_k, cname.str());

    // 2. EDGE INEQUALITIES
    for(int e=0; e < instance->num_edges; e++)
    {
        int u = instance->s[e];
        int v = instance->t[e];

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

    for(int i = 0; i < instance->num_vertices; i++)
        objective_expression += 1*x[i];

    model->setObjective(objective_expression, GRB_MINIMIZE);
    model->update();
}


void Model::create_variables_cc()
{
    char buffer[50];

    // binary selection of each vertex
    x = new GRBVar[instance->num_edges];
    for(int i=0; i < instance->num_edges; i++)
    {
        sprintf(buffer, "x_%d", i);
        x[i] = model->addVar(0.0, 1.0, 1.0, GRB_BINARY, buffer);
    }
    model->update();
}

void Model::create_constraints_cc()
{
    ostringstream cname;

    // 1. FIXED CARDINALITY CONSTRAINT
    GRBLinExpr fixed_cardinality = 0;
    for(int i=0; i < instance->num_edges; i++)
        fixed_cardinality += x[i];

    cname.str("");
    cname << "C1_fixed_cardinality";
    model->addConstr(fixed_cardinality == card_k, cname.str());

    // 2. EDGE INEQUALITIES
    for(int e=0; e < instance->num_conflicts; e++)
    {
        int u = instance->conflicts[e].first;
        int v = instance->conflicts[e].second;

        GRBLinExpr edge_ineq = 0;
        edge_ineq += x[u];
        edge_ineq += x[v];

        cname.str("");
        cname << "C2_edge_" << u << "_" << v;
        model->addConstr(edge_ineq <= 1, cname.str());
    }

    model->update();
}

void Model::create_objective_cc()
{
    GRBLinExpr objective_expression = 0;

    for(int i = 0; i < instance->num_edges; i++)
        objective_expression += (instance->w[i])*x[i];

    model->setObjective(objective_expression, GRB_MINIMIZE);
    model->update();
}

Model::Model(IO *instance, int cardinality)
{
    this->instance = instance;
    this->card_k = cardinality;

    try
    {
        this->env = new GRBEnv();
        this->model = new GRBModel(*env);

        cout << "creating variables..." << endl;
        //create_variables();
        create_variables_cc();

        cout << "creating constraints..." << endl;
        //create_constraints();
        create_constraints_cc();

        cout << "setting objective function..." << endl;
        //create_objective();
        create_objective_cc();

        //model->write("k-stab.lp");
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

int Model::solve()
{
    try
    {
        // turn off all gurobi cut generators
        model->getEnv().set(GRB_IntParam_Cuts, 0);

        // turn off gurobi presolve and heuristics
        model->getEnv().set(GRB_IntParam_Presolve, 0);
        model->getEnv().set(GRB_IntParam_PrePasses, 0);
        model->getEnv().set(GRB_DoubleParam_Heuristics, 0);

        // need to disable presolve reductions that affect user cuts
        model->getEnv().set(GRB_IntParam_PreCrush, 1);

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

int Model::solve_lp_relax()
{
    try
    {
        // turn off all gurobi cut generators
        model->getEnv().set(GRB_IntParam_Cuts, 0);

        // turn off gurobi presolve and heuristics
        model->getEnv().set(GRB_IntParam_Presolve, 0);
        model->getEnv().set(GRB_IntParam_PrePasses, 0);
        model->getEnv().set(GRB_DoubleParam_Heuristics, 0);

        // further options to get a clean LP relaxation
        model->getEnv().set(GRB_IntParam_Method, 0);
        model->getEnv().set(GRB_IntParam_Symmetry, 0);
        model->getEnv().set(GRB_DoubleParam_PreSOS1BigM, 0);
        //model->getEnv().set(GRB_StringParam_ResultFile, "solution.sol");

        //for (int i=0; i < instance->num_vertices; ++i)
        for (int i=0; i < instance->num_edges; ++i)
            x[i].set(GRB_CharAttr_VType, GRB_CONTINUOUS);

        //model->write("mstcc.lp");
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

bool Model::check_half_integer_solution()
{
    if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        for (int i=0; i < instance->num_vertices; ++i)
        {
            double var = x[i].get(GRB_DoubleAttr_X);
            if (var != 0  && var != 0.5  && var != 1)
            {
                //cout << endl << "FOUND A SOLUTION WHICH IS NOT HALF-INTEGER! x[" << i << "] = " << var  << endl;
                return false;
            }
            //cout << x[i].get(GRB_DoubleAttr_X) << endl;
        }
    }

    return true;
}