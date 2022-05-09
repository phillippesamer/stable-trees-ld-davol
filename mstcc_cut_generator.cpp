#include "mstcc_cut_generator.h"

/// algorithm setup switches

bool USE_FAST_FOLKLORE_CUT_IN_IP = false;
bool USE_FAST_INTEGER_CUT_IN_IP = true;

int  SEC_STRATEGY = ALL_CUTS;
bool STORE_SEC_CUT_POOL = false;

#define SEC_ORTHOGONALITY_TOL 0.1
#define SEC_VIOLATION_TOL_IN_IP 0.0001

// at most 14 without changing everything to long double (which gurobi ignores)
#define SEC_SEPARATION_PRECISION_IN_IP 14

SSTCutGenerator::SSTCutGenerator(GRBModel *model, GRBVar *x_vars, IO *instance)
: KStabCutGenerator(model, x_vars, instance)
{
    this->sec_counter = 0;
    this->sec_stats = new sec_statistics();
}

SSTCutGenerator::~SSTCutGenerator()
{
    delete sec_stats;
}

void SSTCutGenerator::callback()
{
    /***
     * The actual callback method within the solver. Currently, only used for 
     * adding cuts/lazy constraints dynamically.
     */

    try
    {
        // callback from the search at a given MIP node: including USER CUTS
        if (where == GRB_CB_MIPNODE)
        {
            // node relaxation solution must be available at the current node
            if (this->getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL)
                return;

            x_val = this->getNodeRel(x_vars, num_vars);

            /***
             * Find violated subtour eliminaton constraints (if any) and/or
             * violated odd-cycle inequalities and add cuts to the model.
             * NB! CURRENTLY SEARCHING FOR OCI'S ONLY WHEN NO VIOLATED SEC
             * WAS EXISTS.
             */
            bool model_updated = run_sec_separation(ADD_USER_CUTS);

            if (!model_updated)
                run_oci_separation(ADD_USER_CUTS);

            delete[] x_val;
        }

        // callback from a new MIP incumbent: including LAZY CONSTRAINTS
        else if (where == GRB_CB_MIPSOL)
        {
            x_val = this->getSolution(x_vars, num_vars);

            // find violated subtour eliminaton constraints (if any) and add cuts to the model
            run_sec_separation(ADD_LAZY_CNTRS);

            delete[] x_val;
        }
    }
    catch (GRBException e)
    {
        cout << "Error " << e.getErrorCode() << " during callback: ";
        cout << e.getMessage() << endl;
    }
    catch (...)
    {
        cout << "Unexpected error during callback" << endl;
    }
}


bool SSTCutGenerator::separate_lpr()
{
    /// Interface to be used when solving the LP relaxation only.

    try
    {
        x_val = new double[num_vars];
        for (long i=0; i < num_vars; ++i)
            x_val[i] = x_vars[i].get(GRB_DoubleAttr_X);

        // find violated subtour eliminaton constraints (if any) and add cuts to the model
        bool model_updated = run_sec_separation(ADD_STD_CNTRS);

        delete[] x_val;
        return model_updated;
    }
    catch (GRBException e)
    {
        cout << "Error " << e.getErrorCode() << " during separate_lpr: ";
        cout << e.getMessage() << endl;
        return false;
    }
    catch (...)
    {
        cout << "Unexpected error during separate_lpr" << endl;
        return false;
    }
}

bool SSTCutGenerator::run_sec_separation(int kind_of_cuts)
{
    /// wrapper for the separation procedure to suit different kinds of cuts

    bool model_updated = false;

    // eventual cuts are stored here
    vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
    vector<long> cuts_rhs = vector<long>();

    if (USE_FAST_INTEGER_CUT_IN_IP)
    {
        // check if solution is integral and attempt more efficient
        // separation procedure in that case
        model_updated = separate_sec_integer(cuts_lhs,cuts_rhs);
    }

    if (!model_updated && USE_FAST_FOLKLORE_CUT_IN_IP)
    {
        // if not using the integer separation or if solution is
        // fractional, try standard separation procedure (capacities 
        // in auxiliary network equal to current vars), which is
        // faster than the classical algorithm but inexact
        model_updated = separate_sec_folklore(cuts_lhs,cuts_rhs);
    }

    if (!model_updated)
    {
        // original 1983 separation procedure (capacities in auxiliary
        // network set so as to yield most violated inequality)
        model_updated = separate_sec_classical(cuts_lhs,cuts_rhs);
    }

    if (model_updated)
    {
        // add cuts
        for (unsigned long idx = 0; idx<cuts_lhs.size(); ++idx)
        {
            if (kind_of_cuts == ADD_USER_CUTS)
                addCut(cuts_lhs[idx] <= cuts_rhs[idx]);

            else if (kind_of_cuts == ADD_LAZY_CNTRS)
                addLazy(cuts_lhs[idx] <= cuts_rhs[idx]);

            else // kind_of_cuts == ADD_STD_CNTRS
                model->addConstr(cuts_lhs[idx] <= cuts_rhs[idx]);
        }
    }

    return model_updated;
}

void SSTCutGenerator::dfs_checking_acyclic( long u,
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

bool SSTCutGenerator::separate_sec_integer( vector<GRBLinExpr> &cuts_lhs,
                                            vector<long> &cuts_rhs )
{
    /***
     * solve the separation problem for an integer solution (in this case, we 
     * may check for violated secs faster)
     */

    // store indices of variables which are set to 0 or 1
    vector<long> s_ = vector<long>();
    for (long i=0; i < num_vars; ++i)
    {
        if (x_val[i] <= EPSILON_TOL || x_val[i] > 1 - EPSILON_TOL)
            s_.push_back(i);
    }

    // proceed only if the current solution is integer
    if (s_.size() == (unsigned) this->num_vars)
    {
        /*
        #ifdef DEBUG_LPR
        cout << "trying fast integer cut... ";
        #endif
        */

        // makes adjlist of the subgraph induced by the current solution
        vector< list<long> > s_adj_list;
        for (long i=0; i<instance->graph->num_vertices; ++i)
            s_adj_list.push_back(list<long>());
        
        for (vector<long>::iterator it = s_.begin(); it != s_.end(); ++it)
        {
            if (x_val[*it] > 1 - EPSILON_TOL)
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
            #ifdef DEBUG_LPR
            cout << "success: ";
            #endif
            */

            long cycle_size = cycle_edges.size();

            // lhs
            while (!cycle_edges.empty())
            {
                violated_constraint += x_vars[cycle_edges.top()];

                /*
                #ifdef DEBUG_LPR
                cout << "+ x[" << cycle_edges.top() << "]";
                #endif
                */

                cycle_edges.pop();
            }

            // save this SEC in the reference arg (added selectively to model by the caller function)
            cuts_lhs.push_back(violated_constraint);
            cuts_rhs.push_back(cycle_size - 1);

            /*
            #ifdef DEBUG_LPR
            cout << " <= " << cycle_size-1 << endl;
            #endif
            */

            return true;
        }

        /*
        #ifdef DEBUG_LPR
        cout << "done" << endl;
        #endif
        */
    }

    return false;
}

bool SSTCutGenerator::separate_sec_folklore( vector<GRBLinExpr> &cuts_lhs,
                                             vector<long> &cuts_rhs )
{
    /***
     * try to solve the separation problem for subtour elimination constraints
     * with a faster (inexact) procedure: just like the regular 1983 method of
     * Padberg and Wolsey, but on a simpler network whose arc capacities are
     * exactly the values of the corresponding vars in the current solution
     */

    // prevent floating point errors by ignoring digits beyond set precision 
    for (long i=0; i < this->num_vars; ++i)
    {
        double tmp = x_val[i] * std::pow(10,SEC_SEPARATION_PRECISION_IN_IP);
        tmp = std::round(tmp);
        x_val[i] = tmp * std::pow(10,-SEC_SEPARATION_PRECISION_IN_IP);
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
    for (long i=0; i<this->num_vars; ++i)
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
                violated_sec *sec = new violated_sec(this->num_vars);
                bool delay_sec_destruction = false;

                GRBLinExpr violated_constraint = 0;

                // lhs: variables corresponding to arcs within checked vertices
                for (long i=0; i<instance->graph->num_vertices; ++i)
                {
                    for (long j=i+1; j<instance->graph->num_vertices; ++j)
                    {
                        long edge_idx = this->instance->graph->index_matrix[i][j];
                        if(edge_idx >= 0 && cutset_chk_v[i] && cutset_chk_v[j])
                        {
                            violated_constraint += x_vars[edge_idx];

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

                        if (counting_cuts.count(sec_id) == 0)
                        {
                            // first time seeing this inequality
                            counting_cuts[sec_id] = 1;

                            if (SEC_STRATEGY == ALL_CUTS)
                            {
                                /*
                                #ifdef DEBUG_LPR
                                cout << "Added cut: " ;
                                for (vector<long>::iterator it = sec->S.begin(); it != sec->S.end(); ++it)
                                    cout << "+ x[" << *it << "] ";
                                cout << "<= " << sec->vertex_count - 1 << endl;
                                #endif
                                */

                                // store sec (caller method adds it to the model)
                                cuts_lhs.push_back(violated_constraint);
                                cuts_rhs.push_back(cutset_size - 1);

                                if (STORE_SEC_CUT_POOL)
                                    sec_stats->pool[sec_id] = 1;
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
                        
                        // lhs: vars corresponding to arcs within unchecked vertices
                        for (long i=0; i<instance->graph->num_vertices; ++i)
                        {
                            for (long j=i+1; j<instance->graph->num_vertices; ++j)
                            {
                                long edge_idx = instance->graph->index_matrix[i][j];
                                if(edge_idx >= 0 && !cutset_chk_v[i] && !cutset_chk_v[j])
                                {
                                    alt_violated_constraint += x_vars[edge_idx];

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

                            if (counting_cuts.count(sec_id) == 0)
                            {
                                // first time seeing this inequality
                                counting_cuts[sec_id] = 1;

                                if (SEC_STRATEGY == ALL_CUTS)
                                {
                                    /*
                                    #ifdef DEBUG_LPR
                                    cout << "Added cut: " ;
                                    for (vector<long>::iterator it = sec->S.begin(); it != sec->S.end(); ++it)
                                        cout << "+ x[" << *it << "] ";
                                    cout << "<= " << sec->vertex_count - 1 << endl;
                                    #endif
                                    */

                                    // store sec (caller method adds it to the model)
                                    cuts_lhs.push_back(alt_violated_constraint);
                                    cuts_rhs.push_back(cutset_size - 1);

                                    if (STORE_SEC_CUT_POOL)
                                        sec_stats->pool[sec_id] = 1;
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
    
    if (separated > 0 && SEC_STRATEGY != ALL_CUTS)
    {
        /*
        #ifdef DEBUG_LPR
        cout << "Added cut(s):" << endl;
        #endif
        */

        // 6.1 ADD MOST VIOLATED CUT, INCLUDED IN BOTH ENHANCED STRATEGIES
        GRBLinExpr violated_constraint = 0;
        
        // lhs: variables corresponding to arcs within this set S
        long edge_count = cuts[most_violated_idx]->S.size();
        for (long e=0; e<edge_count; ++e)
        {
            long idx = cuts[most_violated_idx]->S.at(e);
            violated_constraint += x_vars[idx];

            /*
            #ifdef DEBUG_LPR
                cout << "+ x[" << idx << "] ";
            #endif
            */
        }
        
        // store sec (caller method adds it to the model)
        cuts_lhs.push_back(violated_constraint);
        cuts_rhs.push_back(cuts[most_violated_idx]->vertex_count - 1);

        if (STORE_SEC_CUT_POOL)
            sec_stats->pool[cuts[most_violated_idx]->toString()] = 1;
        
        /*
        #ifdef DEBUG_LPR
        cout << "<= " << cuts[most_violated_idx]->vertex_count - 1 << endl;
        #endif
        */

        // we are done if the strategy is to add just the most violated cut
        
        if (SEC_STRATEGY == ORTHOGONAL_CUTS)
        {
            /***
             * 6.2 ADD ANY OTHER CUT WHOSE HYPERPLANE HAS INNER PRODUCT WITH
             * THAT OF THE MOST VIOLATED CUT CLOSE TO ZERO
             */
            
            // 2-norm of the vector corresponding to most violated inequality
            double norm_v1 = 0.;
            for(long i=0; i < this->num_vars; ++i)
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
                    for(long i=0; i < this->num_vars; ++i)
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
                    for(long i=0; i < this->num_vars; ++i)
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
                    if (inner_prod <= SEC_ORTHOGONALITY_TOL)
                    {
                        GRBLinExpr constr = 0;
                        
                        // lhs: variables corresponding to this cut
                        long e_count = cuts[cut_idx]->S.size();
                        for (long e=0; e<e_count; ++e)
                        {
                            long edge_idx = cuts[cut_idx]->S.at(e);
                            constr += x_vars[edge_idx];
                        }

                        // store sec (caller method adds it to the model)
                        cuts_lhs.push_back(constr);
                        cuts_rhs.push_back(cuts[cut_idx]->vertex_count - 1);

                        if (STORE_SEC_CUT_POOL)
                            sec_stats->pool[cuts[cut_idx]->toString()] = 1;

                        /*
                        #ifdef DEBUG_LPR
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
    for (vector<violated_sec*>::iterator it = cuts.begin(); it != cuts.end(); ++it)
        delete *it;
    cuts.clear();

    return (separated > 0);
}

bool SSTCutGenerator::separate_sec_classical( vector<GRBLinExpr> &cuts_lhs,
                                              vector<long> &cuts_rhs )
{
    /***
     * Solve the separation problem for subtour elimination constraints.
     * Original reference: 1983 Trees and cuts [Padberg, Wolsey]
     */

    // prevent floating point errors by ignoring digits beyond set precision 
    vector<double> x_singleton_cap(instance->graph->num_vertices, 0.);
    for (long i=0; i < this->num_vars; ++i)
    {
        double tmp = x_val[i] * std::pow(10,SEC_SEPARATION_PRECISION_IN_IP);
        tmp = std::round(tmp);
        x_val[i] = tmp * std::pow(10,-SEC_SEPARATION_PRECISION_IN_IP);

        // capacities of singleton cutsets
        long v1 = instance->graph->s[i];
        long v2 = instance->graph->t[i];
        x_singleton_cap[v1] = x_singleton_cap[v1] + x_val[i];
        x_singleton_cap[v2] = x_singleton_cap[v2] + x_val[i];
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
    for (long i=0; i < this->num_vars; ++i)
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

    /***
     * 4. CHECK CUTS IN ORDER OF CAPACITY, ADDING VIOLATED INEQUALITIES (IF ANY)
     * ACCORDING TO CUT INCLUSION STRATEGIES: MOST VIOLATED CUT ONLY, ALL CUTS,
     * OR THE MOST VIOLATED AND THOSE CLOSE ENOUGH TO BEING ORTHOGONAL TO IT
     */

    sort(cuts_values.begin(), cuts_values.end(), sort_pairs_by_snd_val);

    long separated = 0;
    map<string,long> counting_cuts;  // store how many different cuts were found

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
        violated_sec *sec = new violated_sec(this->num_vars);
        bool delay_sec_destruction = false;

        GRBLinExpr violated_constraint = 0;

        // lhs: variables corresponding to edges within checked vertices
        for (long i=0; i<instance->graph->num_vertices; ++i)
        {
            for (long j=i+1; j<instance->graph->num_vertices; ++j)
            {
                long edge_idx = this->instance->graph->index_matrix[i][j];
                if(edge_idx >= 0 && chk_v[i] && chk_v[j])
                {
                    violated_constraint += x_vars[edge_idx];

                    sec->S.push_back(edge_idx);
                    sec->coefficients[edge_idx] = true;

                    // current value of vars, to compute violation
                    cutset_var_sum += x_val[edge_idx];
                }
            }
        }
        
        // 5. VERTICES IN THE CUTSET YIELD VIOLATED SEC
        // TODO: NO NEED FOR A TOLERANCE HERE? if ( cutset_var_sum > ((double) cutset_size - 1.0 + SEC_VIOLATION_TOL_IN_IP) )
        if ( cutset_var_sum > ((double) cutset_size - 1.0) )
        {
            ++separated;  // flag a cut was found

            sec->vertex_count = cutset_size;
            sec->infeasibility = (double) cutset_var_sum 
                    - ((double)cutset_size-1.);

            string sec_id = sec->toString();

            if (counting_cuts.count(sec_id) == 0)
            {
                // first time seeing this inequality
                counting_cuts[sec_id] = 1;

                if (SEC_STRATEGY == ALL_CUTS)
                {
                    /*
                    #ifdef DEBUG_LPR
                    cout << "Added cut: " ;
                    for (vector<long>::iterator it = sec->S.begin(); it != sec->S.end(); ++it)
                        cout << "+ x[" << *it << "] ";
                    cout << "<= " << sec->vertex_count - 1 << endl;
                    #endif
                    */

                    // store sec (caller method adds it to the model)
                    cuts_lhs.push_back(violated_constraint);
                    cuts_rhs.push_back(cutset_size - 1);

                    if (STORE_SEC_CUT_POOL)
                        sec_stats->pool[sec_id] = 1;
                }
                else   // will inspect and add SECs selectively 
                {
                    delay_sec_destruction = true;
                    cuts.push_back(sec);
                }
            }
            else
                counting_cuts[sec_id] = counting_cuts[sec_id] + 1;
        }
        else
        {
            // inequalities are sorted in violation order, so we may stop as
            // soon as a non-violated one is found
            done = true;
        }

        if (!delay_sec_destruction)
            delete sec;

        ++cut_it;

    } // iteration processing each min cut

    /***
     * 6. adds the most violated inequality, as well as those whose hyperplane
     * is "sufficiently close" to being orthogonal with the most violated one
     */

    if (separated > 0 && SEC_STRATEGY != ALL_CUTS)
    {
        // the most violated inequality is given by the first cut above
        long most_violated_idx = 0;

        // 6.1 ADD MOST VIOLATED CUT, INCLUDED IN BOTH ENHANCED STRATEGIES
        GRBLinExpr violated_constraint = 0;
        
        // lhs: variables corresponding to edges within this set S
        long edge_count = cuts[most_violated_idx]->S.size();
        for (long e=0; e<edge_count; ++e)
        {
            long idx = cuts[most_violated_idx]->S.at(e);
            violated_constraint += x_vars[idx];
        }
        
        // store sec (caller method adds it to the model)
        cuts_lhs.push_back(violated_constraint);
        cuts_rhs.push_back(cuts[most_violated_idx]->vertex_count - 1);
        
        if (STORE_SEC_CUT_POOL)
            sec_stats->pool[cuts[most_violated_idx]->toString()] = 1;

        // 6.2 WE ARE DONE IF THE STRATEGY IS TO ADD JUST THE MOST VIOLATED CUT

        if (SEC_STRATEGY == ORTHOGONAL_CUTS)
        {
            // 6.3 ADD ANY OTHER CUT WHOSE HYPERPLANE HAS INNER PRODUCT WITH THAT
            // OF THE MOST VIOLATED CUT CLOSE TO ZERO

            // 2-norm of the vector corresponding to most violated inequality
            double norm_v1 = 0.;
            for(long i=0; i < this->num_vars; ++i)
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
                    for(long i=0; i < this->num_vars; ++i)
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
                    for(long i=0; i < this->num_vars; ++i)
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
                    if (inner_prod <= SEC_ORTHOGONALITY_TOL)
                    {
                        GRBLinExpr constr = 0;
                        
                        // lhs: variables corresponding to this cut
                        long e_count = cuts[cut_idx]->S.size();
                        for (long e=0; e<e_count; ++e)
                        {
                            long edge_idx = cuts[cut_idx]->S.at(e);
                            constr += x_vars[edge_idx];
                        }

                        // store sec (caller method adds it to the model)
                        cuts_lhs.push_back(constr);
                        cuts_rhs.push_back(cuts[cut_idx]->vertex_count - 1);
                        
                        if (STORE_SEC_CUT_POOL)
                            sec_stats->pool[cuts[cut_idx]->toString()] = 1;
                    }
                }

            } // checked each candidate cut
        }

    } // enhanced cut addition strategies

    // clean up
    for (vector<violated_sec*>::iterator it = cuts.begin(); it != cuts.end(); ++it)
        delete *it;
    cuts.clear();

    return (separated > 0);
}
