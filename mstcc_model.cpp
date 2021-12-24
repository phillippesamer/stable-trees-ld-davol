#include "mstcc_model.h"

// cut selection strategies
#define MOST_VIOLATED_CUT 1
#define ORTHOGONAL_CUTS 2
#define ALL_CUTS 3

#define ORTHOGONALITY_TOL 0.1

// algorithm setup switches
int USE_SEC_STRATEGY = ALL_CUTS;
bool USE_INCUMBENT_CHECK = true;
bool USE_FAST_INTEGER_CUT = true;
bool STORE_CUT_POOL = true;

bool sort_pairs_by_snd_val(pair<long,double> a, pair<long,double> b) { return ( a.second<b.second ); }

/// information of a violated subtour elimination constraint
class violated_sec
{
public:
    violated_sec(long edge_count)
    {
        this->edge_count = edge_count;
        coefficients = vector<bool>(edge_count, false);
    }
    virtual ~violated_sec()
    {
        //
    }
    
    string toString()
    {
        ostringstream sec_lhs;
        sec_lhs.str("");
        for (vector<bool>::iterator it = coefficients.begin(); it != coefficients.end(); ++it)
            sec_lhs << *it;
        return sec_lhs.str();
    }

     void reset()
     {
        S.clear();
        vertex_count = 0;
        infeasibility = 0.;
        coefficients = vector<bool>(edge_count, false);
     }

    long edge_count;            // instance parameter
    vector<long> S;             // index of edge vars in cut
    long vertex_count;          // number of vertices spanned by S
    double infeasibility;       // amount by which the sec is violated
    vector<bool> coefficients;  // hyperplane coefficients
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
            #endif

            // eventual cuts are stored here
            vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
            vector<long> cuts_rhs = vector<long>();

            if (USE_FAST_INTEGER_CUT)
            {
                // more efficient separation procedure for an integral solution
                model_updated = separate_SEC_integer(cuts_lhs,cuts_rhs);
            }

            /*
            if (!model_updated)
            {
                // if not using the integer separation or if solution is
                // fractional, try standard separation procedure
                model_updated = separate_SEC(cuts_lhs,cuts_rhs);
            }
            */

            if (!model_updated)
            {
                // TODO: replace call above by this corrected separation procedure
                model_updated = separate_SEC_fallback(cuts_lhs,cuts_rhs);
            }
            

            if (model_updated)
            {
                //cout << cuts_lhs.size() << " cuts added" << endl;

                // add cut(s)
                for (unsigned long idx = 0; idx<cuts_lhs.size(); ++idx)
                    model->addConstr(cuts_lhs[idx] <= cuts_rhs[idx]);

                // reoptmize
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

        /*
        #ifdef DEBUG
        if (STORE_CUT_POOL)
            cout << stats->pool.size() << " SEC cuts added" << endl;

        cout << "fractional vars in the final solution of the lpr";
        bool any_frac_var = false;
        for (long i=0; i < this->instance->graph->num_edges; ++i)
        {
            double tmp = this->x[i].get(GRB_DoubleAttr_X);
            if (tmp > EPSILON_TOL && tmp < 1 - EPSILON_TOL)
            {
                if (!any_frac_var)
                {
                    cout << endl;
                    any_frac_var = true;
                }

                cout.precision(10);
                cout << "\tx[ (" << instance->graph->s[i] << "," << instance->graph->t[i] << ") ] = " << this->x[i].get(GRB_DoubleAttr_X) << endl;
            }
        }
        if (!any_frac_var)
            cout << "... none!" << endl;

        cout << "runtime (s): " << this->lp_runtime << endl;
        #endif
        */

// TODO: REMOVE THIS
instance->graph->mst();
cout << instance->instance_id << "\t"
     << instance->graph->mst_weight - model->get(GRB_DoubleAttr_ObjVal) << "\t"
     << this->lp_runtime << endl;
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

        double tmp2 = tmp1 * std::pow(10,SEPARATION_PRECISION);
        tmp2 = std::round(tmp2);
        tmp2 = tmp2 * std::pow(10,-SEPARATION_PRECISION);

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
    
    SmartDigraph lemon_g;
    SmartDigraph::ArcMap<double> lemon_cap(lemon_g);

    vector<SmartDigraph::Node> lemon_vertices;
    for (long i=0; i<instance->graph->num_vertices; ++i)
        lemon_vertices.push_back(lemon_g.addNode());

    vector< vector<SmartDigraph::Arc> > lemon_arcs;
    for (long i=0; i<instance->graph->num_vertices; ++i)
        lemon_arcs.push_back( vector<SmartDigraph::Arc>(instance->graph->num_vertices, lemon::Invalid()) );
    
    // for each edge in the original graph, create arcs on both directions
    for (long i=0; i<instance->graph->num_edges; ++i)
    {
        long v1 = instance->graph->s[i];
        long v2 = instance->graph->t[i];
        
        // LEMON arcs
        lemon_arcs[v1][v2] = lemon_g.addArc(lemon_vertices[v1], lemon_vertices[v2]);
        lemon_cap[ lemon_arcs[v1][v2] ] = x_val[i];

        lemon_arcs[v2][v1] = lemon_g.addArc(lemon_vertices[v2], lemon_vertices[v1]);
        lemon_cap[lemon_arcs[v2][v1]] = x_val[i];
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

    for (long sink=0; sink<instance->graph->num_vertices; ++sink)
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
                vector<bool> cutset_chk_v = vector<bool>(instance->graph->num_vertices, false);
                
                for (long idx=0; idx<instance->graph->num_vertices; ++idx)
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
                bool delay_sec_destruction = false;

                GRBLinExpr violated_constraint = 0;

                // lhs: variables (x) corresponding to arcs within checked vertices
                for (long i=0; i<instance->graph->num_vertices; ++i)
                {
                    for (long j=i+1; j<instance->graph->num_vertices; ++j)
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
                
                // 5.1 VERTICES IN THE CUTSET YIELD VIOLATED SEC
                // ONLY ACCEPTING A CUT IF VIOLATION IS BEYOND A TOLERANCE HERE
                if ( cutset_var_sum > ((double) cutset_size - 1.0 + EPSILON_TOL) )
                {
                    if (cutset_size == 0 || cutset_size == instance->graph->num_vertices)
                        cerr << endl << endl << "\t\tUNEXPECTED ERROR: NON PROPER SUBSET S IN A SEC" << endl;
                    else
                    {
                        ++separated;  // flag a cut was found

                        sec->vertex_count = cutset_size;
                        sec->infeasibility = (double) cutset_var_sum 
                                - ((double)cutset_size-1.);

                        // counting different cuts (only relevant when including all/orthogonal cuts)
                        string sec_id = sec->toString();

                        if (counting_cuts[sec_id] <= 0)
                        {
                            // first time seeing this inequality
                            counting_cuts[sec_id] = 1;


                            if (USE_SEC_STRATEGY == ALL_CUTS)
                            {
                                /*
                                #ifdef DEBUG
                                cout << "Added cut: " ;
                                for (vector<long>::iterator it = sec->S.begin(); it != sec->S.end(); ++it)
                                    cout << "+ x[" << *it << "] ";
                                cout << "<= " << sec->vertex_count - 1 << endl;
                                #endif
                                */

                                // store sec (caller method adds it to the model)
                                cuts_lhs.push_back(violated_constraint);
                                cuts_rhs.push_back(cutset_size - 1);

                                if (STORE_CUT_POOL)
                                    stats->pool[sec_id] = 1;
                            }
                            else   // will inspect/add secs selectively
                            {
                                delay_sec_destruction = true;
                                cuts.push_back(sec);

                                if (sec->infeasibility > most_violated)
                                {
                                    most_violated = sec->infeasibility;
                                    most_violated_idx = cuts.size() - 1;
                                }
                            }
                        }
                        else
                            counting_cuts[sec_id] = counting_cuts[sec_id] + 1;
                    }
                }

                // 5.2 VERTICES IN THE COMPLEMENT OF THE CUTSET YIELDS VIOLATED SEC
                else
                {
                    cutset_size = instance->graph->num_vertices - cutset_size;
                    if (cutset_size == 0 || cutset_size == instance->graph->num_vertices)
                        cerr << endl << endl << "\t\tUNEXPECTED ERROR: NON PROPER SUBSET S\' IN A SEC" << endl;
                    else
                    {
                        sec->reset();

                        GRBLinExpr alt_violated_constraint = 0;
                        cutset_var_sum = 0;
                        
                        // lhs: vars (x) corresponding to arcs within unchecked vertices
                        for (long i=0; i<instance->graph->num_vertices; ++i)
                        {
                            for (long j=i+1; j<instance->graph->num_vertices; ++j)
                            {
                                long edge_idx = instance->graph->index_matrix[i][j];
                                if(edge_idx >= 0 && !cutset_chk_v[i] && !cutset_chk_v[j])
                                {
                                    alt_violated_constraint += x[edge_idx];

                                    sec->S.push_back(edge_idx);
                                    sec->coefficients[edge_idx] = true;
                                    
                                    // current value of vars, to compute violation
                                    cutset_var_sum += x_val[edge_idx];
                                }
                            }
                        }
                        
                        // ONLY ACCEPTING A CUT IF VIOLATION IS BEYOND A TOLERANCE HERE
                        if ( cutset_var_sum > ((double) cutset_size - 1.0 + EPSILON_TOL) )
                        {
                            ++separated;  // flag a cut was found

                            sec->vertex_count = cutset_size;
                            sec->infeasibility = (double) cutset_var_sum 
                                - ((double)cutset_size-1.);

                            // counting different cuts (only relevant when including all/orthogonal cuts)
                            string sec_id = sec->toString();

                            if (counting_cuts[sec_id] <= 0)
                            {
                                // first time seeing this inequality
                                counting_cuts[sec_id] = 1;

                                if (USE_SEC_STRATEGY == ALL_CUTS)
                                {
                                    /*
                                    #ifdef DEBUG
                                    cout << "Added cut: " ;
                                    for (vector<long>::iterator it = sec->S.begin(); it != sec->S.end(); ++it)
                                        cout << "+ x[" << *it << "] ";
                                    cout << "<= " << sec->vertex_count - 1 << endl;
                                    #endif
                                    */

                                    // store sec (caller method adds it to the model)
                                    cuts_lhs.push_back(alt_violated_constraint);
                                    cuts_rhs.push_back(cutset_size - 1);

                                    if (STORE_CUT_POOL)
                                        stats->pool[sec_id] = 1;
                                }
                                else   // will inspect/add secs selectively
                                {
                                    delay_sec_destruction = true;
                                    cuts.push_back(sec);

                                    if (sec->infeasibility > most_violated)
                                    {
                                        most_violated = sec->infeasibility;
                                        most_violated_idx = cuts.size() - 1;
                                    }
                                }
                            }
                            else
                                counting_cuts[sec_id] = counting_cuts[sec_id] + 1;
                        }
                        //else
                        //    cerr << endl << endl << "\t\t WARNING: MINCUT < 1 BUT COULD NOT FIND VIOLATED SEC! (LIKELY DUE TO VIOLATION TOLERANCE)" << endl;
                    }
                }

                if (!delay_sec_destruction)
                    delete sec;

            }   // done with a single mincut subproblem
                
        }

    } // |V|-1 min cut computations
    
    // 6. ENHANCED CUT ADDITION STRATEGIES
    
    if (separated > 0 && USE_SEC_STRATEGY != ALL_CUTS)
    {
        /*
        #ifdef DEBUG
        cout << "Added cut(s):" << endl;
        #endif
        */

        // counting different cuts
        stats->sec_diff_cuts.push_back(counting_cuts.size());
        
        // 6.1 ADD MOST VIOLATED CUT, INCLUDED IN BOTH ENHANCED STRATEGIES
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

        if (STORE_CUT_POOL)
            stats->pool[cuts[most_violated_idx]->toString()] = 1;
        
        /*
        #ifdef DEBUG
        cout << "<= " << cuts[most_violated_idx]->vertex_count - 1 << endl;
        #endif
        */

        // we are done if the strategy is to add just the most violated cut
        
        if (USE_SEC_STRATEGY == ORTHOGONAL_CUTS)
        {
            // 6.2 ADD ANY OTHER CUT WHOSE HYPERPLANE HAS INNER PRODUCT WITH THAT
            // OF THE MOST VIOLATED CUT CLOSE TO ZERO
            
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
                if (cut_idx != (unsigned) most_violated_idx)
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

                        if (STORE_CUT_POOL)
                            stats->pool[cuts[cut_idx]->toString()] = 1;

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
    //for(vector<violated_sec*>::iterator it=cuts.begin(); it!=cuts.end(); ++it)
    //    delete (*it);
    cuts.clear();

    return (separated > 0);
}


bool StableSpanningTreeModel::separate_SEC_fallback( vector<GRBLinExpr> &cuts_lhs,
                                                     vector<long> &cuts_rhs )
{
    /***
     * Solve the separation problem for subtour elimination constraints.
     * Original reference: 1983 Trees and cuts [Padberg, Wolsey]
     */

    vector<double> x_val;
    vector<double> x_singleton_cap(instance->graph->num_vertices, 0.);
    for (long i=0; i < instance->graph->num_edges; ++i)
    {
        // current solution values, up to the precision set
        /*
        if (tmp < EPSILON_TOL)
            x_val.push_back(0.0);
        else if (tmp > 1.0 - EPSILON_TOL)
            x_val.push_back(1.0);
        else
            x_val.push_back(tmp);
        */
        double tmp1 = this->x[i].get(GRB_DoubleAttr_X);

        double tmp2 = tmp1 * std::pow(10,SEPARATION_PRECISION);
        tmp2 = std::round(tmp2);
        tmp2 = tmp2 * std::pow(10,-SEPARATION_PRECISION);

        x_val.push_back(tmp2);

        // capacities of singleton cutsets
        long v1 = instance->graph->s[i];
        long v2 = instance->graph->t[i];
        x_singleton_cap[v1] = x_singleton_cap[v1] + tmp2;
        x_singleton_cap[v2] = x_singleton_cap[v2] + tmp2;
    }

    // LEMON digraph representing the current solution
    
    SmartDigraph lemon_g;
    SmartDigraph::ArcMap<double> lemon_cap(lemon_g);

    vector<SmartDigraph::Node> lemon_vertices;
    for (long i=0; i<instance->graph->num_vertices; ++i)
        lemon_vertices.push_back(lemon_g.addNode());

    vector< vector<SmartDigraph::Arc> > lemon_arcs;
    for (long i=0; i<instance->graph->num_vertices; ++i)
        lemon_arcs.push_back( vector<SmartDigraph::Arc>(instance->graph->num_vertices, lemon::Invalid()) );

    // for each edge in the original graph, create arcs on both directions
    for (long i=0; i<instance->graph->num_edges; ++i)
    {
        long v1 = instance->graph->s[i];
        long v2 = instance->graph->t[i];
        
        // LEMON arcs
        if (x_val[i] > 0)
        {
            lemon_arcs[v1][v2] = lemon_g.addArc(lemon_vertices[v1], lemon_vertices[v2]);
            lemon_cap[ lemon_arcs[v1][v2] ] = 0.5 * x_val[i];

            lemon_arcs[v2][v1] = lemon_g.addArc(lemon_vertices[v2], lemon_vertices[v1]);
            lemon_cap[lemon_arcs[v2][v1]] = 0.5 * x_val[i];
        }
    }

    // arcs from dummy vertex s to each other vertex in the original graph
    SmartDigraph::Node lemon_dummy_s = lemon_g.addNode();
    vector<SmartDigraph::Arc> lemon_arcs_from_s(instance->graph->num_vertices, lemon::Invalid());
    for (long i=0; i<instance->graph->num_vertices; ++i)
    {
        lemon_arcs_from_s[i] = lemon_g.addArc(lemon_dummy_s, lemon_vertices[i]);
        lemon_cap[ lemon_arcs_from_s[i] ] = std::max( 0.5 * x_singleton_cap[i] - 1 , 0.0 );
    }

    // arcs from each vertex in the original graph to dummy vertex t
    SmartDigraph::Node lemon_dummy_t = lemon_g.addNode();
    vector<SmartDigraph::Arc> lemon_arcs_to_t(instance->graph->num_vertices, lemon::Invalid());
    for (long i=0; i<instance->graph->num_vertices; ++i)
    {
        lemon_arcs_to_t[i] = lemon_g.addArc(lemon_vertices[i], lemon_dummy_t);
        lemon_cap[ lemon_arcs_to_t[i] ] = std::max( 1 - 0.5 * x_singleton_cap[i] , 0.0 );
    }
    
    /***
     * 1. |V| max flow (min cut) computations: from the dummy source to the
     * dummy sink, each forcing k \in V to lie in the cutset by making the
     * capacity of arc (s,k) arbitrarily large
     */
    
    vector<violated_sec*> cuts;
    vector< pair<long,double> > cuts_values;   // <fixed vertex, cut capacity>
    vector<long> cuts_cutset_sizes;
    vector< vector<bool> > cuts_cutsets;

    for (long fixed_vertex=0; fixed_vertex<instance->graph->num_vertices; ++fixed_vertex)
    {
        // update capacity of the arc (source,fixed) to "infinity"
        // N.B. |V|^2 i a simple upper bound to the sum of all capacities
        lemon_cap[ lemon_arcs_from_s[fixed_vertex] ] = 
            std::pow(instance->graph->num_vertices, 2); //numeric_limits<int>::max();

        /***
         * 2. LEMON: compute maxflow/mincut value using first phase of
         * Goldberg & Tarjan preflow push-relabel algorithm (with "highest 
         * label" and "bound decrease" heuristics). The worst case time 
         * complexity of the algorithm is in O(n^2 * m^0.5)
         */

        Preflow<SmartDigraph, SmartDigraph::ArcMap<double> > lemon_preflow(
            lemon_g, lemon_cap, lemon_dummy_s, lemon_dummy_t);

        lemon_preflow.runMinCut();
        double cut_val = lemon_preflow.flowValue();

        // 3. RETRIEVE AND STORE CUTSET S
        long cutset_size = 0;
        vector<bool> chk_v = vector<bool>(instance->graph->num_vertices, false);

        for (long idx=0; idx<instance->graph->num_vertices; ++idx)
        {
            // query if node is on the source side of the minimum cut found
            if ( lemon_preflow.minCut(lemon_vertices[idx]) )
            {
                chk_v[idx] = true;
                ++cutset_size;
            }
        }

        cuts_values.push_back( make_pair(fixed_vertex,cut_val) );
        cuts_cutset_sizes.push_back(cutset_size);
        cuts_cutsets.push_back(chk_v);

        // restore capacity of the arc (source,fixed)
        lemon_cap[ lemon_arcs_from_s[fixed_vertex] ] 
            = std::max( 0.5 * x_singleton_cap[fixed_vertex] - 1 , 0.0 );

    } // |V| min cut computations

    // inspect cutsets in order of increasing capacity
    sort(cuts_values.begin(), cuts_values.end(), sort_pairs_by_snd_val);

    /***
     * 4. CHECK CUTS IN ORDER OF CAPACITY, ADDING VIOLATED INEQUALITIES (IF ANY)
     * ACCORDING TO CUT INCLUSION STRATEGIES: MOST VIOLATED CUT ONLY, ALL CUTS,
     * OR THE MOST VIOLATED AND THOSE CLOSE ENOUGH TO BEING ORTHOGONAL TO IT
     */

    double most_violated = -1;
    long most_violated_idx = -1;
    map<string,long> counting_cuts;  // store how many different cuts were found

    long separated = 0;
    bool done = false;
    vector<pair<long,double> >::iterator cut_it = cuts_values.begin();

    while ( cut_it != cuts_values.end() && !done )
    {
        // for each cut, inspect whether the cutset yields a violated SEC 

        long fixed_vertex = cut_it->first;
        long cutset_size = cuts_cutset_sizes[fixed_vertex];
        vector<bool>& chk_v = cuts_cutsets[fixed_vertex];

        double cutset_var_sum = 0;
        
        // object storing new cut
        violated_sec *sec = new violated_sec(instance->graph->num_edges);
        bool delay_sec_destruction = false;

        GRBLinExpr violated_constraint = 0;

        // lhs: variables (x) corresponding to edges within checked vertices
        for (long i=0; i<instance->graph->num_vertices; ++i)
        {
            for (long j=i+1; j<instance->graph->num_vertices; ++j)
            {
                long edge_idx = this->instance->graph->index_matrix[i][j];
                if(edge_idx >= 0 && chk_v[i] && chk_v[j])
                {
                    violated_constraint += x[edge_idx];

                    sec->S.push_back(edge_idx);
                    sec->coefficients[edge_idx] = true;

                    // current value of vars, to compute violation
                    cutset_var_sum += x_val[edge_idx];
                }
            }
        }
        
        // 5.1 VERTICES IN THE CUTSET YIELD VIOLATED SEC
        // TODO: NO NEED FOR A TOLERANCE HERE? if ( cutset_var_sum > ((double) cutset_size - 1.0 + VIOLATION_TOL) )
        if ( cutset_var_sum > ((double) cutset_size - 1.0) )
        {
            if (cutset_size == 0 || cutset_size == instance->graph->num_vertices)
                cerr << endl << endl << "\t\tUNEXPECTED ERROR: NON PROPER SUBSET S IN A SEC" << endl;
            else
            {
                ++separated;  // flag a cut was found

                sec->vertex_count = cutset_size;
                sec->infeasibility = (double) cutset_var_sum 
                        - ((double)cutset_size-1.);

                // counting different cuts (only relevant when including all/orthogonal cuts)
                string sec_id = sec->toString();

                if (counting_cuts[sec_id] <= 0)
                {
                    // first time seeing this inequality
                    counting_cuts[sec_id] = 1;

                    if (USE_SEC_STRATEGY == ALL_CUTS)
                    {
                        /*
                        #ifdef DEBUG
                        cout << "Added cut: " ;
                        for (vector<long>::iterator it = sec->S.begin(); it != sec->S.end(); ++it)
                            cout << "+ x[" << *it << "] ";
                        cout << "<= " << sec->vertex_count - 1 << endl;
                        #endif
                        */

                        // store sec (caller method adds it to the model)
                        cuts_lhs.push_back(violated_constraint);
                        cuts_rhs.push_back(cutset_size - 1);

                        if (STORE_CUT_POOL)
                            stats->pool[sec_id] = 1;
                    }
                    else   // will inspect and add SECs selectively 
                    {
                        delay_sec_destruction = true;
                        cuts.push_back(sec);

                        if (sec->infeasibility > most_violated)
                        {
                            most_violated = sec->infeasibility;
                            most_violated_idx = cuts.size() - 1;
                        }
                    }
                }
                else
                    counting_cuts[sec_id] = counting_cuts[sec_id] + 1;
            }
        }

        // 5.2 VERTICES IN THE COMPLEMENT OF THE CUTSET YIELDS VIOLATED SEC
        else
        {
            done = true;
        }
        /*
        {
            cutset_size = instance->graph->num_vertices - cutset_size;
            if (cutset_size == 0 || cutset_size == instance->graph->num_vertices)
                cerr << endl << endl << "\t\tUNEXPECTED ERROR: NON PROPER SUBSET S\' IN A SEC" << endl;
            else
            {
                sec->reset();

                GRBLinExpr alt_violated_constraint = 0;
                cutset_var_sum = 0;
                
                // lhs: vars (x) corresponding to edges within unchecked vertices
                for (long i=0; i<instance->graph->num_vertices; ++i)
                {
                    for (long j=i+1; j<instance->graph->num_vertices; ++j)
                    {
                        long edge_idx = instance->graph->index_matrix[i][j];
                        if(edge_idx >= 0 && !chk_v[i] && !chk_v[j])
                        {
                            alt_violated_constraint += x[edge_idx];

                            sec->S.push_back(edge_idx);
                            sec->coefficients[edge_idx] = true;
                            
                            // current value of vars, to compute violation
                            cutset_var_sum += x_val[edge_idx];
                        }
                    }
                }
                
                // TODO: NO NEED FOR A TOLERANCE HERE? if ( cutset_var_sum > ((double) cutset_size - 1.0 + VIOLATION_TOL) )
                if ( cutset_var_sum > ((double) cutset_size - 1.0) )
                {
cout.precision(12);
cout << endl << endl << "HAD TO INCLUDE SEC FROM COMPLEMENT" << cutset_var_sum << " > |\\bar S|-1 = " << (double) cutset_size - 1.0 << endl;

                    ++separated;  // flag a cut was found

                    sec->vertex_count = cutset_size;
                    sec->infeasibility = (double) cutset_var_sum 
                        - ((double)cutset_size-1.);

cout << "violation: " << sec->infeasibility << endl;

                    // counting different cuts (only relevant when including all/orthogonal cuts)
                    string sec_id = sec->toString();

                    if (counting_cuts[sec_id] <= 0)
                    {
                        // first time seeing this inequality
                        counting_cuts[sec_id] = 1;

                        if (USE_SEC_STRATEGY == ALL_CUTS)
                        {
                            //#ifdef DEBUG
                            //cout << "Added cut: " ;
                            //for (vector<long>::iterator it = sec->S.begin(); it != sec->S.end(); ++it)
                            //    cout << "+ x[" << *it << "] ";
                            //cout << "<= " << sec->vertex_count - 1 << endl;
                            //#endif

                            // store sec (caller method adds it to the model)
                            cuts_lhs.push_back(alt_violated_constraint);
                            cuts_rhs.push_back(cutset_size - 1);

                            if (STORE_CUT_POOL)
                                stats->pool[sec_id] = 1;
                        }
                        else   // will inspect and add SECs selectively
                        {
                            delay_sec_destruction = true;
                            cuts.push_back(sec);

                            if (sec->infeasibility > most_violated)
                            {
                                most_violated = sec->infeasibility;
                                most_violated_idx = cuts.size() - 1;
                            }
                        }
                    }
                    else
                        counting_cuts[sec_id] = counting_cuts[sec_id] + 1;
                }
            }
        }
        */

        if (!delay_sec_destruction)
            delete sec;

        ++cut_it;

    } // iteration processing each min cut

    /***
     * 6. adds the most violated inequality, as well as those whose hyperplane
     * is "sufficiently close" to being orthogonal with the most violated one
     */

    if (separated > 0 && USE_SEC_STRATEGY != ALL_CUTS)
    {
        /*
        #ifdef DEBUG
        cout << "Added cut(s):" << endl;
        #endif
        */

        // counting different cuts
        stats->sec_diff_cuts.push_back(counting_cuts.size());
        
        // 6.1 ADD MOST VIOLATED CUT, INCLUDED IN BOTH ENHANCED STRATEGIES
        GRBLinExpr violated_constraint = 0;
        
        // lhs: variables (x) corresponding to edges within this set S
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
        
        if (STORE_CUT_POOL)
            stats->pool[cuts[most_violated_idx]->toString()] = 1;

        /*
        #ifdef DEBUG
        cout << "<= " << cuts[most_violated_idx]->vertex_count - 1 << endl;
        #endif
        */

        // we are done if the strategy is to add just the most violated cut

        if (USE_SEC_STRATEGY == ORTHOGONAL_CUTS)
        {
            // 6.2 ADD ANY OTHER CUT WHOSE HYPERPLANE HAS INNER PRODUCT WITH THAT
            // OF THE MOST VIOLATED CUT CLOSE TO ZERO

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
                if (cut_idx != (unsigned) most_violated_idx)
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
                        
                        if (STORE_CUT_POOL)
                            stats->pool[cuts[cut_idx]->toString()] = 1;

                        /*
                        #ifdef DEBUG
                        for (vector<long>::iterator it = cuts[cut_idx]->S.begin(); it != cuts[cut_idx]->S.end(); ++it)
                            cout << "+ x[" << *it << "] ";
                        cout << "<= " << cuts[cut_idx]->vertex_count - 1 << endl;
                        #endif
                        */
                    }
                }

            } // checked each candidate cut
        }

    } // enhanced cut addition strategies

    // clean up
    //for(vector<violated_sec*>::iterator it=cuts.begin(); it!=cuts.end(); ++it)
    //    delete (*it);
    cuts.clear();

    return (separated > 0);
}
