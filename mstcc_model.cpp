#include "mstcc_model.h"

// cut selection strategies
#define FIRST_CUT 0
#define MOST_VIOLATED_CUT 1
#define ORTHOGONAL_CUTS 2
#define ALL_CUTS 3

#define ORTHOGONALITY_TOL 0.1

// algorithm setup switches
int USE_SEC_STRATEGY = ALL_CUTS;
bool USE_INCUMBENT_CHECK = true;
bool USE_FAST_INTEGER_CUT = true;

/// information of a violated subtour elimination constraint
class violated_sec
{
public:
    violated_sec(long edge_count)
    {
        this->edge_count = edge_count;
        coefficients = new bool[edge_count];
        memset(coefficients, 0, sizeof(bool)*edge_count);
    }
    virtual ~violated_sec()
    {
        delete[] coefficients;
    }
    
    string toString()
    {
        long len = edge_count+1;
        char buffer[len];
        memset(buffer, 0, sizeof(char)*len);
        
        for (long i=0; i<len; ++i)
        {
            if (i == len-1)
                buffer[i] = '\0';
            else
                buffer[i] = (coefficients[i] == 0) ? '0' : '1';
        }
        
        string id(buffer);
        return id;
    }
    
    // index all vars (to make it easier to check orthogonality)
    bool *coefficients;
    
    long edge_count;           // instance parameter
    vector<long> S;            // index of edge vars in cut
    long vertex_count;         // number of vertices spanned by S
    double infeasibility;      // amount by which the sec is violated
};


StableSpanningTreeModel::StableSpanningTreeModel(IO *instance)
: KStabModel(instance)
{
    this->lp_bound = this->lp_runtime = -1;
    stats = new statistics();

    // add redundant constraints (all vertices should have degree >= 1)
    // to make LP relaxation faster
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


StableSpanningTreeModel::~StableSpanningTreeModel()
{
    delete stats;
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
        /*
        model->set(GRB_IntParam_Presolve, 0);
        model->set(GRB_DoubleParam_PreSOS1BigM, 0);
        model->set(GRB_DoubleParam_PreSOS2BigM, 0);
        model->set(GRB_IntParam_PreSparsify, 0);
        model->set(GRB_IntParam_PreCrush, 1);
        model->set(GRB_IntParam_DualReductions, 0);
        model->set(GRB_IntParam_Aggregate, 0);
        */

        if (logging == true)
            model->set(GRB_IntParam_OutputFlag, 1);
        else
            model->set(GRB_IntParam_OutputFlag, 0);

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
            model_updated = false;

            #ifdef DEBUG
            cout << "LP relaxation pass #" << lp_passes << " (bound = "
                 << model->get(GRB_DoubleAttr_ObjVal) << ")" << endl;
            
            /*
            for (long i=0; i < this->instance->graph->num_edges; ++i)
                cout << "x[" << i << "] = " << this->x[i].get(GRB_DoubleAttr_X) << endl;
            */

            
            cout << "fractional vars in the final solution of the lpr" << endl;
            for (long i=0; i < this->instance->graph->num_edges; ++i)
            {
                double tmp = this->x[i].get(GRB_DoubleAttr_X);
                if (tmp > EPSILON_TOL && tmp < 1 - EPSILON_TOL)
                {
                    cout.precision(10);
                    cout << "\tx[ (" << instance->graph->s[i] << "," << instance->graph->t[i] << ") ] = " << this->x[i].get(GRB_DoubleAttr_X) << endl;
                }
            }
            
            #endif

            // eventual cuts are stored here
            vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
            vector<long> cuts_rhs = vector<long>();

            if (USE_FAST_INTEGER_CUT)
            {
                // more efficient separation procedure for an integral solution
                model_updated = separate_SEC_integer(cuts_lhs,cuts_rhs);
            }

            if (!model_updated)
            {
                // if not using the integer separation or if solution is
                // fractional, try standard separation procedure
                model_updated = separate_SEC(cuts_lhs,cuts_rhs);
            }

/*            if (!model_updated)
            {
                // if the most efficient separation procedures failed, we make a
                // last attempt checking the subgraph induced by fractional vars
                model_updated = separate_SEC_fallback(cuts_lhs,cuts_rhs);
            }
*/
            if (model_updated)
            {
                // add cut(s)
                for (unsigned long idx = 0; idx<cuts_lhs.size(); ++idx)
                {
                    // TODO: filter cuts here?
                    model->addConstr(cuts_lhs[idx] <= cuts_rhs[idx]);
                }

                // reoptmize
                model->update();
                model->optimize();
                this->lp_passes++;
                cout << "model status = " << model->get(GRB_IntAttr_Status) << endl;
            }
        }

        // LP relaxation time
        gettimeofday(clock_stop, 0);
        unsigned long clock_time = 1.e6 * (clock_stop->tv_sec - clock_start->tv_sec) +
                                          (clock_stop->tv_usec - clock_start->tv_usec);
        this->lp_runtime = ((double)clock_time / (double)1.e6);
        free(clock_start);
        free(clock_stop);

// TODO: REMOVE THIS
instance->graph->mst();
cout << instance->graph->mst_weight << endl;
cout << model->get(GRB_DoubleAttr_ObjVal) << endl;
cout << instance->graph->mst_weight - model->get(GRB_DoubleAttr_ObjVal) << endl;
//model->write("mstcc_with_sec.lp");

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


void StableSpanningTreeModel::dfs_checking_acyclic( long u,
                                                    long parent,
                                                    vector<bool> &check,
                                                    long &checked_count,
                                                    stack<long> &cycle_edges,
                                                    vector< list<long> > &adj_list,
                                                    bool &acyclic)
{
    /// specialized depth-first search to check if a subgraph is acyclic

    check[u] = true;
    ++checked_count;

    list<long>::iterator it = adj_list[u].begin();
    while (it != adj_list[u].end() && acyclic)
    {
        long v = (*it);
        
        long edge_idx = this->instance->graph->index_matrix[u][v];

        if (!check[v])
        {
            cycle_edges.push(edge_idx);
            dfs_checking_acyclic(v, u, check, checked_count, cycle_edges, adj_list, acyclic);
            if (acyclic)
                cycle_edges.pop();
        }
        else if (v != parent)
        {
            // vertex (not the antecessor) already visited: back edge!
            acyclic = false;
            cycle_edges.push(edge_idx);
        }

        ++it;
    }
}


bool StableSpanningTreeModel::separate_SEC_integer( vector<GRBLinExpr> &cuts_lhs,
                                                    vector<long> &cuts_rhs )
{
    /***
     * solve the separation problem for an integer solution (in this case, we 
     * may check for violated secs in reduced complexity)
     */

    // store indices of variables which are set to 0 or 1
    vector<long> s_ = vector<long>();
    for (long i=0; i < this->instance->graph->num_edges; ++i)
    {
        double x_val = this->x[i].get(GRB_DoubleAttr_X);
        if (x_val <= EPSILON_TOL || x_val > 1 - EPSILON_TOL)
            s_.push_back(i);
    }

    // proceed only if the current solution is integer
    if (s_.size() == (unsigned) instance->graph->num_edges)
    {
        /*
        #ifdef DEBUG
        cout << "trying fast integer cut... ";
        #endif
        */

        // makes adjlist of the subgraph induced by the current solution
        vector< list<long> > s_adj_list;
        for (long i=0; i<instance->graph->num_vertices; ++i)
            s_adj_list.push_back(list<long>());
        
        for (vector<long>::iterator it = s_.begin(); it != s_.end(); ++it)
        {
            if (this->x[*it].get(GRB_DoubleAttr_X) > 1 - EPSILON_TOL)
            {
                long u = this->instance->graph->s[*it];
                long v = this->instance->graph->t[*it];
                
                s_adj_list[u].push_back(v);
                s_adj_list[v].push_back(u);
            }
        }

        bool acyclic = true;

        stack<long> cycle_edges;
        vector<bool> check = vector<bool>(instance->graph->num_vertices, false);
        long checked_count = 0;

        // the search starts only once, if the subgraph is connected/acyclic
        long root = 0;
        while (root < instance->graph->num_vertices && acyclic && checked_count < instance->graph->num_vertices)
        {
            if (!check[root])
                dfs_checking_acyclic(root, -1, check, checked_count, cycle_edges, s_adj_list, acyclic);

            ++root;
        }

        if (!acyclic)
        {
            // set S with checked vertices yields violated constraint
            GRBLinExpr violated_constraint = 0;
            
            /*
            #ifdef DEBUG
            cout << "success: ";
            #endif
            */

            long cycle_size = cycle_edges.size();

            // lhs
            while (!cycle_edges.empty())
            {
                violated_constraint += x[cycle_edges.top()];

                /*
                #ifdef DEBUG
                cout << "+ x[" << cycle_edges.top() << "]";
                #endif
                */

                cycle_edges.pop();
            }

            // save this SEC in the reference arg (added selectively to model by the caller function)
            cuts_lhs.push_back(violated_constraint);
            cuts_rhs.push_back(cycle_size - 1);

            /*
            #ifdef DEBUG
            cout << " <= " << cycle_size-1 << endl;
            #endif
            */

            return true;
        }

        /*
        #ifdef DEBUG
        cout << "done" << endl;
        #endif
        */
    }

    return false;
}


bool StableSpanningTreeModel::separate_SEC( vector<GRBLinExpr> &cuts_lhs,
                                            vector<long> &cuts_rhs )
{
    /***
     * solve the separation problem for subtour elimination constraints
     */

    vector<double> x_val;
    for (long i=0; i < instance->graph->num_edges; ++i)
    {
        double tmp1 = this->x[i].get(GRB_DoubleAttr_X);

        double tmp2 = tmp1 * 1e10;
        tmp2 = std::round(tmp2);
        tmp2 = tmp2 * 1e-10;

        x_val.push_back(tmp2);

        /*
        if (tmp < EPSILON_TOL)
            x_val.push_back(0.0);
        else if (tmp > 1.0 - EPSILON_TOL)
            x_val.push_back(1.0);
        else
            x_val.push_back(tmp);
        */
    }

    // LEMON digraph representing the current solution
    long num_vertices = instance->graph->num_vertices;
    
    SmartDigraph lemon_g;
    SmartDigraph::ArcMap<double> lemon_cap(lemon_g);

    vector<SmartDigraph::Node> lemon_vertices;
    for (long i=0; i<instance->graph->num_vertices; ++i)
        lemon_vertices.push_back(lemon_g.addNode());
    
    // for each edge in the original graph, create arcs on both directions
    for (long i=0; i<instance->graph->num_edges; ++i)
    {
        long v1 = instance->graph->s[i];
        long v2 = instance->graph->t[i];
        
        // LEMON arcs
        SmartDigraph::Arc lemon_arc_1;
        lemon_arc_1 = lemon_g.addArc(lemon_vertices[v1], lemon_vertices[v2]);
        lemon_cap[lemon_arc_1] = x_val[i];
        
        SmartDigraph::Arc lemon_arc_2;
        lemon_arc_2 = lemon_g.addArc(lemon_vertices[v2], lemon_vertices[v1]);
        lemon_cap[lemon_arc_2] = x_val[i];
    }
    
    // |V|-1 max flow (min cut) computations: from source vertex 0 to k \in V\{source}
    long source = 0;
    
    long separated = 0;
    double most_violated = -1;
    long most_violated_idx = -1;
    vector<violated_sec*> cuts;
    
    // store how many different cuts were generated in this execution
    map<string,long> counting_cuts;
    
    vector<double> mincut_vals;

    for (long sink=0; sink<num_vertices; ++sink)
    {
        if (sink != source)
        {
            /***
             * 1. LEMON: compute maxflow/mincut value using first phase of
             * Goldberg & Tarjan preflow push-relabel algorithm (with "highest 
             * label" and "bound decrease" heuristics). The worst case time 
             * complexity of the algorithm is in O(n^2 * m^0.5)
             */

            Preflow<SmartDigraph, SmartDigraph::ArcMap<double> > lemon_preflow(
                lemon_g, lemon_cap, lemon_vertices[source], lemon_vertices[sink]);

            lemon_preflow.runMinCut();
            
            double mincut = lemon_preflow.flowValue();
            
            mincut_vals.push_back(mincut);

            // 2. THERE EXISTS A VIOLATED SEC IFF A MAX FLOW (MIN CUT) < 1
            if (mincut < 1)
            {
                // 3. LEMON: A MIN CUT IS ALREADY COMPUTED .: RETRIEVE CUTSET S
                long cutset_size = 0;
                vector<bool> cutset_chk_v = vector<bool>(num_vertices, false);
                
                for (long idx=0; idx<num_vertices; ++idx)
                {
                    // query if node is on the source side of the minimum cut found
                    if ( lemon_preflow.minCut(lemon_vertices[idx]) )
                    {
                        cutset_chk_v[idx] = true;
                        ++cutset_size;
                    }
                }
                
                // 4. FIND VIOLATED SEC (DECIDE WHETHER TO USE S OR ITS COMPLEMENT)
                double cutset_var_sum = 0;
                
                // object storing new cut
                violated_sec *sec = new violated_sec(instance->graph->num_edges);

                GRBLinExpr violated_constraint = 0;

                // lhs: variables (x) corresponding to arcs within checked vertices
                for (long i=0; i<num_vertices; ++i)
                {
                    for (long j=i+1; j<num_vertices; ++j)
                    {
                        long edge_idx = this->instance->graph->index_matrix[i][j];
                        if(edge_idx >= 0 && cutset_chk_v[i] && cutset_chk_v[j])
                        {
                            violated_constraint += x[edge_idx];

                            sec->S.push_back(edge_idx);
                            sec->coefficients[edge_idx] = true;

                            // current value of vars, to compute violation
                            cutset_var_sum += x_val[edge_idx];
                        }
                    }
                }
                
                // checked vertices yield sec
                if ( cutset_var_sum > ((double) cutset_size - 1.0 - VIOLATION_TOL) )
                {
                    if (cutset_size == 0 || cutset_size == num_vertices)
                        cerr << endl << endl << "\t\tUNEXPECTED ERROR: NON PROPER SUBSET S IN A SEC" << endl;
                    else
                    {
                        ++separated;  // flag a cut was found

                        sec->vertex_count = cutset_size;
                        sec->infeasibility = (double) cutset_var_sum 
                                - ((double)cutset_size-1.);

                        if (USE_SEC_STRATEGY == FIRST_CUT)
                        {
                            /*
                            #ifdef DEBUG
                            cout << "Added cut: " ;
                            for (vector<long>::iterator it = sec->S.begin(); it != sec->S.end(); ++it)
                                cout << "+ x[" << *it << "] ";
                            cout << "<= " << sec->vertex_count - 1 << endl;
                            #endif
                            */

                            delete sec;

                            // store sec (caller method adds it to the model)
                            cuts_lhs.push_back(violated_constraint);
                            cuts_rhs.push_back(cutset_size - 1);
                        }
                        else   // will compute all secs
                        {
                            cuts.push_back(sec);

                            if (sec->infeasibility > most_violated)
                            {
                                most_violated = sec->infeasibility;
                                most_violated_idx = cuts.size() - 1;
                            }
                            
                            // counting different cuts
                            string sec_id = sec->toString();
                            counting_cuts[sec_id] = counting_cuts[sec_id] <= 0 ?
                                1 : counting_cuts[sec_id] + 1;
                        }
                    }
                }
                else   // unchecked vertices (complement of S) yield sec
                {
                    cutset_size = num_vertices - cutset_size;
                    if (cutset_size == 0 || cutset_size == num_vertices)
                        cerr << endl << endl << "\t\tUNEXPECTED ERROR: NON PROPER SUBSET S\' IN A SEC" << endl;
                    else
                    {
                        delete sec;
                        
                        violated_sec *sec2 = new violated_sec(instance->graph->num_edges);
                        GRBLinExpr alt_violated_constraint = 0;
                        cutset_var_sum = 0;
                        
                        // lhs: vars (x) corresponding to arcs within unchecked vertices
                        for (long i=0; i<num_vertices; ++i)
                        {
                            for (long j=i+1; j<num_vertices; ++j)
                            {
                                long edge_idx = instance->graph->index_matrix[i][j];
                                if(edge_idx >= 0 && !cutset_chk_v[i] && !cutset_chk_v[j])
                                {
                                    alt_violated_constraint += x[edge_idx];

                                    sec2->S.push_back(edge_idx);
                                    sec2->coefficients[edge_idx] = true;
                                    
                                    // current value of vars, to compute violation
                                    cutset_var_sum += x_val[edge_idx];
                                }
                            }
                        }
                        
                        if ( cutset_var_sum > ((double) cutset_size - 1.0 - VIOLATION_TOL) )
                        {
                            ++separated;  // flag a cut was found

                            sec2->vertex_count = cutset_size;
                            sec2->infeasibility = (double) cutset_var_sum 
                                - ((double)cutset_size-1.);

                            if (USE_SEC_STRATEGY == FIRST_CUT)
                            {
                                /*
                                #ifdef DEBUG
                                cout << "Added cut: " ;
                                for (vector<long>::iterator it = sec2->S.begin(); it != sec2->S.end(); ++it)
                                    cout << "+ x[" << *it << "] ";
                                cout << "<= " << sec2->vertex_count - 1 << endl;
                                #endif
                                */

                                delete sec2;

                                // store sec (caller method adds it to the model)
                                cuts_lhs.push_back(alt_violated_constraint);
                                cuts_rhs.push_back(cutset_size - 1);
                            }
                            else   // will compute all secs
                            {
                                cuts.push_back(sec2);

                                if (sec2->infeasibility > most_violated)
                                {
                                    most_violated = sec2->infeasibility;
                                    most_violated_idx = cuts.size() - 1;
                                }
                                
                                // counting different cuts
                                string sec_id = sec2->toString();
                                counting_cuts[sec_id] = counting_cuts[sec_id] <= 0 ? 
                                    1 : counting_cuts[sec_id] + 1;
                            }
                        }
                        else
                            cerr << endl << endl << "\t\tUNEXPECTED ERROR: MINCUT < 1 BUT COULD NOT FIND VIOLATED SEC!" << endl;
                    }
                }

            }   // if a cut was found in this single mincut subproblem, we stored it
                
            // 5. FINISHING ONE MAXFLOW COMPUTATION
            
            // return after finding first violated sec?
            if (separated > 0 && USE_SEC_STRATEGY == FIRST_CUT)
            {
                vector<violated_sec*>::iterator it;                                                                                 
                for(it=cuts.begin(); it!=cuts.end(); ++it)
                    delete (*it);
                cuts.clear();
                
                return true;
            }
        }

    } // |V|-1 min cut computations
    
    // 6. ENHANCED CUT ADDITION STRATEGIES
    
    long sec_count = 0;
    
    if (separated > 0 && USE_SEC_STRATEGY != FIRST_CUT)
    {
        /*
        #ifdef DEBUG
        cout << "Added cut(s):" << endl;
        #endif
        */

        // counting different cuts
        stats->sec_diff_cuts.push_back(counting_cuts.size());
        
        // 6.1 ADD MOST VIOLATED CUT, INCLUDED IN ALL FOLLOWING STRATEGIES
        GRBLinExpr violated_constraint = 0;
        
        // lhs: variables (x) corresponding to arcs within this set S
        long edge_count = cuts[most_violated_idx]->S.size();
        for (long e=0; e<edge_count; ++e)
        {
            long idx = cuts[most_violated_idx]->S.at(e);
            violated_constraint += x[idx];

            /*
            #ifdef DEBUG
                cout << "+ x[" << idx << "] ";
            #endif
            */
        }
        
        // store sec (caller method adds it to the model)
        cuts_lhs.push_back(violated_constraint);
        cuts_rhs.push_back(cuts[most_violated_idx]->vertex_count - 1);

        /*
        #ifdef DEBUG
        cout << "<= " << cuts[most_violated_idx]->vertex_count - 1 << endl;
        #endif
        */

        ++sec_count;
        
        // store cut identifier (directing vector) to avoid adding same cuts
        // (except for the most violated one, which we add without checking!)
        stats->pool[cuts[most_violated_idx]->toString()] = 1;
        
        // 6.2 WE ARE DONE IF (USE_SEC_STRATEGY == MOST_VIOLATED_CUT)
        
        if (USE_SEC_STRATEGY == ALL_CUTS)
        {
            for (unsigned cut_idx=0; cut_idx<cuts.size(); ++cut_idx)
            {
                if ( cut_idx != (unsigned) most_violated_idx  &&
                     stats->pool[cuts[cut_idx]->toString()] != 1 )  // not a repeated cut
                {
                    GRBLinExpr constr = 0;
                    
                    // lhs: variables (x) corresponding to this cut
                    long e_count = cuts[cut_idx]->S.size();
                    for (long e=0; e<e_count; ++e)
                    {
                        long edge_idx = cuts[cut_idx]->S.at(e);
                        constr += x[edge_idx];
                    }
                    
                    // store sec (caller method adds it to the model)
                    cuts_lhs.push_back(constr);
                    cuts_rhs.push_back(cuts[cut_idx]->vertex_count - 1);

                    // store cut identifier (directing vector) to avoid adding same cut
                    stats->pool[cuts[cut_idx]->toString()] = 1;

                    ++sec_count;

                    /*
                    #ifdef DEBUG
                    for (vector<long>::iterator it = cuts[cut_idx]->S.begin(); it != cuts[cut_idx]->S.end(); ++it)
                        cout << "+ x[" << *it << "] ";
                    cout << "<= " << cuts[cut_idx]->vertex_count - 1 << endl;
                    #endif
                    */
                }
            }
        }
        else if (USE_SEC_STRATEGY == ORTHOGONAL_CUTS)
        {
            // ADD EVERY OTHER CUT WHICH IS "SUFFICIENTLY ORTHOGONAL" TO THE FIRST
            
            // 2-norm of the vector corresponding to most violated inequality
            double norm_v1 = 0.;
            for(long i=0; i<instance->graph->num_edges; ++i)
            {
                double base = (double) cuts[most_violated_idx]->coefficients[i];
                double sq = pow(base, 2.);
                norm_v1 += sq;
            }
            //double rhs1 = cuts[most_violated_idx]->vertex_count - 1;
            //norm_v1 += pow(rhs1, 2);   // rhs: |S| - 1
            norm_v1 = sqrt(norm_v1);
            
            for (unsigned cut_idx=0; cut_idx<cuts.size(); ++cut_idx)
            {
                if ( cut_idx != (unsigned) most_violated_idx &&
                     stats->pool[cuts[cut_idx]->toString()] != 1 )  // not a repeated cut
                {
                    // 2-norm of candidate cut vector
                    double norm_v2 = 0.;
                    for(long i=0; i<instance->graph->num_edges; ++i)
                    {
                        double base = (double) cuts[cut_idx]->coefficients[i];
                        double sq = pow(base, 2.);
                        norm_v2 += sq;
                    }
                    //double rhs2 = cuts[cut_idx]->vertex_count - 1;
                    //norm_v2 += pow(rhs2, 2);   // rhs: |S| - 1
                    norm_v2 = sqrt(norm_v2);
                
                    // inner product
                    double dot = 0.;
                    for(long i=0; i<instance->graph->num_edges; ++i)
                    {
                        double v1_i = (double) cuts[most_violated_idx]->coefficients[i];
                        double v2_i = (double) cuts[cut_idx]->coefficients[i];
                        dot += (v1_i * v2_i);
                    }
                    //dot += (rhs1*rhs2);  // adds because both rhs are non-negative
                    
                    // add cut if normalized product is close to 0
                    double norm = norm_v1 * norm_v2;
                    double inner_prod = (double) dot / norm;
                    
                    // add sec if sufficiently close to being orthogonal to the most violated one
                    if (inner_prod <= ORTHOGONALITY_TOL)
                    {
                        GRBLinExpr constr = 0;
                        
                        // lhs: variables (x) corresponding to this cut
                        long e_count = cuts[cut_idx]->S.size();
                        for (long e=0; e<e_count; ++e)
                        {
                            long edge_idx = cuts[cut_idx]->S.at(e);
                            constr += x[edge_idx];
                        }

                        // store sec (caller method adds it to the model)
                        cuts_lhs.push_back(constr);
                        cuts_rhs.push_back(cuts[cut_idx]->vertex_count - 1);
                        
                        // store cut identifier (directing vector) to avoid adding same cut
                        stats->pool[cuts[cut_idx]->toString()] = 1;

                        ++sec_count;

                        /*
                        #ifdef DEBUG
                        for (vector<long>::iterator it = cuts[cut_idx]->S.begin(); it != cuts[cut_idx]->S.end(); ++it)
                            cout << "+ x[" << *it << "] ";
                        cout << "<= " << cuts[cut_idx]->vertex_count - 1 << endl;
                        #endif
                        */
                    }
                }

            } // close to orthogonal candidate cuts stored

        }

    } // enhanced cut addition strategies

    // clean up
    vector<violated_sec*>::iterator it;
    for(it=cuts.begin(); it!=cuts.end(); ++it)
        delete (*it);
    cuts.clear();

    for (vector<double>::iterator it = mincut_vals.begin(); it < mincut_vals.end(); ++it)
        cout << *it << " ";
    cout << endl;

    return (sec_count > 0);
}


bool StableSpanningTreeModel::separate_SEC_fallback( vector<GRBLinExpr> &cuts_lhs,
                                                     vector<long> &cuts_rhs )
{
    /***
     * try to separate solution with the cutset of vertices induced by
     * edges with fractional values
     */

    vector<double> x_val;
    for (long i=0; i < instance->graph->num_edges; ++i)
    {
        double tmp = this->x[i].get(GRB_DoubleAttr_X);
        if (tmp < EPSILON_TOL)
            x_val.push_back(0.0);
        else if (tmp > 1.0 - EPSILON_TOL)
            x_val.push_back(1.0);
        else
            x_val.push_back(tmp);
    }

    set<long> frac_vars;
    cout << "fractional vars in the final solution of the lpr" << endl;
    for (long i=0; i < this->instance->graph->num_edges; ++i)
    {
        double tmp = this->x[i].get(GRB_DoubleAttr_X);
        if (tmp > EPSILON_TOL && tmp < 1 - EPSILON_TOL)
        {
            frac_vars.insert(instance->graph->s[i]);
            frac_vars.insert(instance->graph->t[i]);

            cout.precision(10);
            cout << "\tx[ (" << instance->graph->s[i] << "," << instance->graph->t[i] << ") ] = " << this->x[i].get(GRB_DoubleAttr_X) << endl;
        }
    }

    cout << "the subgraph they induce:" << endl;
    for (set<long>::iterator it = frac_vars.begin(); it != frac_vars.end(); ++it)
    {
        cout << "\tN(" << *it << ") = { ";
        list<long>::iterator neighbour = instance->graph->adj_list[*it].begin();
        while (neighbour != instance->graph->adj_list[*it].end())
        {
            cout << *neighbour << " ";
            long edge_idx = instance->graph->index_matrix[*it][*neighbour];
            if (edge_idx>0)
                cout << "(" <<  setw(3) << this->x[edge_idx].get(GRB_DoubleAttr_X) << "   ) \t";
            else
                cout << "(   -   ) \t";
            ++neighbour;
        }
        cout << "}" << endl;
    }

    long cutset_size = frac_vars.size();

    vector<bool> cutset_chk_v = vector<bool>(instance->graph->num_vertices, false);
    for (set<long>::iterator it = frac_vars.begin(); it != frac_vars.end(); ++it)
        cutset_chk_v[*it] = true;
    
    double cutset_var_sum = 0;
    
    vector<long> cutset_induced_edges;

    GRBLinExpr violated_constraint = 0;

    // lhs: variables (x) corresponding to edges within checked vertices
    for (long i=0; i<instance->graph->num_vertices; ++i)
    {
        for (long j=i+1; j<instance->graph->num_vertices; ++j)
        {
            long edge_idx = this->instance->graph->index_matrix[i][j];
            if(edge_idx >= 0 && cutset_chk_v[i] && cutset_chk_v[j])
            {
                violated_constraint += x[edge_idx];

                cutset_induced_edges.push_back(edge_idx);

                // current value of vars, to compute violation
                cutset_var_sum += x_val[edge_idx];
            }
        }
    }
    
    cout << "sum = " << cutset_var_sum << " (while |S|-1= " << cutset_size - 1 << ")" << endl;

    // checked vertices yield sec
    if ( cutset_var_sum > ((double) cutset_size - 1.0 - VIOLATION_TOL) )
    {
        if (cutset_size == 0 || cutset_size == instance->graph->num_vertices)
            cerr << endl << endl << "\t\tUNEXPECTED ERROR: NON PROPER SUBSET S IN A SEC" << endl;
        else
        {
            #ifdef DEBUG
            cout << "Added cut: " ;
            for (vector<long>::iterator it = cutset_induced_edges.begin(); it != cutset_induced_edges.end(); ++it)
                cout << "+ x[" << *it << "] ";
            cout << "<= " << cutset_size - 1 << endl;
            #endif

            // store sec (caller method adds it to the model)
            cuts_lhs.push_back(violated_constraint);
            cuts_rhs.push_back(cutset_size - 1);
            return true;
        }
    }

    return false;
}
