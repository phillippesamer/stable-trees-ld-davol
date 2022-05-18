#include "kstab_model.h"

// algorithm setup switch
bool MAXIMAL_CLIQUE_ENUMERATION = true;

KStabModel::KStabModel(IO *instance)
{
    this->instance = instance;
    this->fixed_cardinality = instance->graph->num_vertices-1;
    this->clique_counter = 0;
    this->clique_sizes = map<long,long>();

    this->solution_weight = numeric_limits<double>::max();
    this->solution_dualbound = numeric_limits<double>::max();
    this->solution_status = STATUS_UNKNOWN;
    this->solution_runtime = -1;

    try
    {
        this->env = new GRBEnv();
        this->model = new GRBModel(*env);

        create_variables();
        create_constraints();
        create_objective();

        cutgen = new KStabCutGenerator(model, x, instance);

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
    delete cutgen;
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

    if (MAXIMAL_CLIQUE_ENUMERATION)
    {
        /***
         * 2. ADDING ALL MAXIMAL CLIQUE INEQUALITIES
         * NB! We enumerate maximal cliques IN THE CONFLICT GRAPH! Below we
         * define num_vertices and num_edges for the sake of clarity only!
         */

        long num_vertices = instance->graph->num_edges;
        long num_edges = instance->num_conflicts;

        // optimized representation of the search tree and adjacencies in the algorithm
        vector<bool> nodes(num_vertices, true);   // "subg" in Tomita et al. (2006)
        vector<bool> open(num_vertices, true);    // "cand" in Tomita et al. (2006)

        // adjacency matrix ("gamma" in Tomita et al. 2006)
        this->cliques_bit_adj = vector<vector<bool> >();
        for (long i=0; i<num_vertices; ++i)
            this->cliques_bit_adj.push_back( vector<bool>(num_vertices,false) );

        for (long c=0; c<num_edges; ++c)
        {
            pair<long,long> endpoints = instance->conflicts.at(c);
            long v1 = endpoints.first;
            long v2 = endpoints.second;
            
            this->cliques_bit_adj[v1][v2] = true;
            this->cliques_bit_adj[v2][v1] = true;
        }

        // 3. run maximal cliques enumeration algorithm (Tomita et al. 2006)
        GRBLinExpr initial_constraint = 0;
        all_maximal_cliques(nodes, open, initial_constraint);

        // done
        this->cliques_bit_adj.clear();

        cout << clique_counter << " clique inequalities added ";
        if (clique_counter < instance->num_conflicts)
            cout << "(LESS than the " << instance->num_conflicts << " edge inequalities)";
        else if (clique_counter > instance->num_conflicts)
            cout << "(MORE than the " << instance->num_conflicts << " edge inequalities)";
        cout << endl;

        for (map<long,long>::iterator it = clique_sizes.begin(); it != clique_sizes.end(); ++it)
        {
            if (it->second == 1)
                cout << setw(6) << it->second << " " << it->first << "-clique" << endl;
            else
                cout << setw(6) << it->second << " " << it->first << "-cliques" << endl;
        }
    }
    else
    {
        // 2. USE SIMPLE EDGE INEQUALITIES
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
    }

    model->update();
}

void KStabModel::all_maximal_cliques(vector<bool> &nodes,
                                     vector<bool> &open,
                                     GRBLinExpr &current_constraint)
{
    /// implements the "cliques" algorithm of Tomita, Tanaka & Takahashi (2006)
    
    // convenience only (enumerating cliques IN THE CONFLICT GRAPH)
    long num_vertices = instance->graph->num_edges;

    // 1. finish a branch in the recursion tree if the current clique is maximal
    bool nodes_allfalse = true;
    vector<bool>::iterator nodes_it = nodes.begin();
    while (nodes_it != nodes.end() && nodes_allfalse)
    {
        if (*nodes_it)
            nodes_allfalse = false;

        ++nodes_it;
    }

    if (nodes_allfalse)
    {
        unsigned len = current_constraint.size();

        // no need to include the 1-cliques
        if (len > 1)
        {
            this->model->addConstr(current_constraint <= 1);
    
            ++clique_counter;
            this->clique_sizes[len] = clique_sizes[len] <= 0 ? 1 : clique_sizes[len] + 1;
        }

        return;
    }
    
    // 2. else: we may extend the current clique. First, choose vertex u with most open neighbors
    long max_val = -1;
    long max_node = -1;

    for (long u=0; u<num_vertices; ++u)
    {
        if (nodes[u])
        {
            // count entries in open \cap adj[u]
            long tmp_count = 0;
            for (long i=0; i<num_vertices; ++i)
                if (open[i] && cliques_bit_adj[u][i])
                    ++tmp_count;

            if (tmp_count > max_val)
            {
                max_val = tmp_count;
                max_node = u;
            }
        }
    }

    // 3. ext_u corresponds to non-neighbors of u which are still open, i.e. open - adj[u]
    //vector<bool> ext_u(num_vertices, false);
    vector<long> ext_u_idx;
    for (long i=0; i<num_vertices; ++i)
    {
        if ( open[i] && !(cliques_bit_adj[max_node][i]) )
        {
            //ext_u[i] = true;
            ext_u_idx.push_back(i);
        }
    }

    // 4. one recursive call for each vertex through which the clique can be extended
    for (vector<long>::iterator it = ext_u_idx.begin(); it != ext_u_idx.end(); ++it)
    {
        // vertex q joins the clique
        long q = *it;
        open[q] = false;

        // copy constraint expression, and extend it with q
        GRBLinExpr constraint_q = GRBLinExpr(current_constraint);
        constraint_q += x[q];

        // new nodes and open sets restricting to neighbors of current vertex q
        vector<bool> nodes_q(num_vertices, false);
        vector<bool> open_q(num_vertices, false);
        for (long i=0; i<num_vertices; ++i)
        {
            if ( nodes[i] && cliques_bit_adj[q][i] )
                nodes_q[i] = true;

            if ( open[i] && cliques_bit_adj[q][i] )
                open_q[i] = true;
        }

        // recursive call from the child node
        all_maximal_cliques(nodes_q, open_q, constraint_q);
    }
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

        // should disable presolve reductions that affect user cuts
        model->set(GRB_IntParam_PreCrush, 1);

        // must set parameter indicating presence of lazy constraints
        //model->set(GRB_IntParam_LazyConstraints, 1);

        // set callback to separate OCIs and solve IP
        model->setCallback(this->cutgen);
        model->optimize();

        // additional information on oci cuts added
        if (logging && cutgen->oci_counter > 0)
        {
            cout << cutgen->oci_counter << " odd-cycle inequalities added" << endl;

            for (map<long,long>::iterator it = cutgen->oci_stats->oci_len.begin();
                 it != cutgen->oci_stats->oci_len.end(); ++it)
            {
                if (it->second == 1)
                    cout << setw(6) << it->second << " " << it->first << "-cycle" << endl;
                else
                    cout << setw(6) << it->second << " " << it->first << "-cycles" << endl;
            }
        }

        //model->write("kstab.lp");
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
        cout << "Unexpected error during optimization inside KStabModel::solve()" << endl;
        return 0;
    }
}

int KStabModel::save_optimization_status()
{
    /// Set class fields accordingly after call to optmize()

    this->solution_runtime = model->get(GRB_DoubleAttr_Runtime);

    // store results in this object
    if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        this->solution_status = AT_OPTIMUM;
        this->solution_vector.clear();

        this->solution_weight = this->solution_dualbound 
                              = model->get(GRB_DoubleAttr_ObjVal);

        /*
        #ifdef DEBUG
            cout << "Model optimal vector: " << endl;
        #endif
        */

        // save bool vector of this solution
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

        this->solution_weight = numeric_limits<double>::max();
        this->solution_dualbound = numeric_limits<double>::max();
        this->solution_vector.clear();

        cout << "Model infeasible! (runtime "
             << solution_runtime << ")" << endl;

        return 0;
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
    {
        this->solution_status = STATUS_UNKNOWN;

        this->solution_weight = (model->get(GRB_IntAttr_SolCount) > 0) ?
                                model->get(GRB_DoubleAttr_ObjVal) :
                                numeric_limits<double>::max();
        this->solution_dualbound = model->get(GRB_DoubleAttr_ObjBound);

        cout << "Time limit exceeded (" << solution_runtime << ")" << endl;
        cout << "Dual bound " << this->solution_dualbound 
             << ", primal bound " << this->solution_weight 
             << " (MIP gap " << 100*model->get(GRB_DoubleAttr_MIPGap) << "%)" 
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

double KStabModel::runtime()
{
    return model->get(GRB_DoubleAttr_Runtime);
}

void KStabModel::set_time_limit(double tl)
{
    model->set(GRB_DoubleParam_TimeLimit, tl);
}

void KStabModel::update_single_weight(long idx, double new_weight)
{
    x[idx].set(GRB_DoubleAttr_Obj, new_weight);
    model->update();

    /*
    #ifdef DEBUG
        model->write("kstab_2.lp");
    #endif
    */
}

void KStabModel::update_all_weights(vector<double> new_weights)
{
    GRBLinExpr objective_expression = 0;

    for(long i = 0; i < instance->graph->num_edges; i++)
        objective_expression += (new_weights[i])*x[i];

    model->setObjective(objective_expression, GRB_MINIMIZE);
    model->update();
}

pair<ModelStatus,double> KStabModel::probe_var(long probe_idx, bool probe_value)
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
    double result = numeric_limits<double>::max();

    if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL)
    {
        status = AT_OPTIMUM;
        result = model->get(GRB_DoubleAttr_ObjVal);
    }
    else if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
    {
        status = IS_INFEASIBLE;
        result = numeric_limits<double>::max();
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
