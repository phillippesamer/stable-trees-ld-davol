#include "kstab_cut_generator.h"

/// algorithm setup switches


bool SEPARATE_OCI = true;

int  OCI_STRATEGY = ORTHOGONAL_CUTS;
bool OCI_STRATEGY_SHIFT_AFTER_ROOT = true;
int  OCI_STRATEGY_SHIFT_CHOICE = MOST_VIOLATED_CUT;

bool STORE_OCI_CUT_POOL = false;

bool CUTS_AT_ROOT_ONLY = false;

#define OCI_ORTHOGONALITY_TOL 0.01
#define OCI_VIOLATION_TOL_IN_IP 0

// at most 14 without changing everything to long double (which gurobi ignores)
#define OCI_SEPARATION_PRECISION_IN_IP 14

///////////////////////////////////////////////////////////////////////////////

/***
 * The following implements a binary min-heap based priority queue, used in
 * Dijkstra's SSSP algorithm, which in turn is used in the separation
 * algorithm for OCIs
 */

typedef struct
{
    double estimate;   // d in Dijkstra algorithm in CLRS (2009)
    long v;            // handle to vertex, to get its adjacency list
    long previous;     // key from vertices in the sp to v, including 'source' and 'v')
    long heap_pos;     // handle to current position in the heap
} sssp_heap_element;

typedef struct
{
public:
    void build_min_heap(vector<sssp_heap_element*> &nodes)
    {
        /// constructs min-heap from a vector

        unsigned long length = nodes.size();
        heap_size = length;

        // fill heap (pointer-) vector
        heap.clear();
        heap.push_back(0);   // dummy head
        for (unsigned long k=0; k<length; ++k)
            heap.push_back(nodes[k]);

        // heap property
        for (unsigned long i = floor(length/2); i > 0; --i)
            min_heapify(i);
    }

    sssp_heap_element* min()
    {
        return heap[1];
    }

    sssp_heap_element* extract_min()
    {
        // when exporting a heap interface, signal error for heap underflow:
        if (heap_size < 1) return 0;

        sssp_heap_element *min = heap[1];

        // maintain heap property
        heap[1] = heap[heap_size];
        heap[1]->heap_pos = 1;
        heap[heap_size] = 0;
        --heap_size;
        min_heapify(1);

        min->heap_pos = 0;   // invalidate handle between heap and graph objects
        return min;
    }

    void insert(sssp_heap_element* e)
    {
        double new_key = e->estimate;
        ++heap_size;
        if (heap.size() <= heap_size)    // no empty node (considering dummy head)
            heap.push_back(0);

        // start at new leaf and move parent down until a smaller parent is found
        unsigned long i = heap_size;
        unsigned long parent = floor(i/2);
        while (i > 1 && heap[parent]->estimate > new_key)
        {
            heap[i] = heap[parent];
            heap[i]->heap_pos = i;
            i = parent;
            parent = floor(i/2);
        }

        heap[i] = e;
        heap[i]->heap_pos = i;
    }

    bool decrease_key(unsigned long i, double new_key)
    {
        /// decreases key and return true, or return false case the new key is greater than current one

        // when exporting a heap interface, signal error here
        if (heap[i]->estimate < new_key) return false;

        // start at current node and move toward the root until a smaller parent is found
        unsigned long parent = floor(i/2);
        while (i > 1 && heap[parent]->estimate > new_key)
        {
            // exchange parent <-> current node
            sssp_heap_element *tmp = heap[i];
            heap[i] = heap[parent];
            heap[i]->heap_pos = i;
            heap[parent] = tmp;
            heap[parent]->heap_pos = parent;

            i = parent;
            parent = floor(i/2);
        }

        heap[i]->estimate = new_key;
        return true;
    }

    unsigned long get_size()
    {
        return heap_size;
    }

private:
    void min_heapify(unsigned long root)
    {
        /// ensure min-heap property (starting from element at 'root')

        unsigned long left = 2*root;       // left child
        unsigned long right = 2*root + 1;  // right child

        // determines the smallest among the root and its children
        unsigned long smallest = root;
        if (right > heap_size)
        {
            if (left > heap_size)
                return;
            else
                smallest = left;
        }
        else
        {
            if (heap[left]->estimate <= heap[right]->estimate)
                smallest = left;
            else
                smallest = right;
        }

        // if heap property violated, reorder nodes and continue above
        if (heap[root]->estimate > heap[smallest]->estimate)
        {
            // exchange root
            sssp_heap_element *tmp = heap[root];
            heap[root] = heap[smallest];
            heap[root]->heap_pos = root;
            heap[smallest] = tmp;
            heap[smallest]->heap_pos = smallest;

            min_heapify(smallest);
        }
    }

    unsigned long heap_size;
    vector<sssp_heap_element*> heap;

} sssp_minheap;

///////////////////////////////////////////////////////////////////////////////

/// specialized depth-first search to identify/count connected components

void dfs_to_tag_components(long u,
                           long component, 
                           vector<long> &chk, 
                           vector< vector<long> > &adj_list)
{
    // auxiliary procedure implementing dfs to check connected components in adj_list
    chk[u] = component;

    for (unsigned i=0; i<adj_list[u].size(); ++i)
    {
        long v = adj_list[u].at(i);
        if (chk[v] < 0)
            dfs_to_tag_components(v, component, chk, adj_list);
    }
}

vector<long> check_components(vector< vector<long> > &adj_list)
{
    long num_vertices = adj_list.size();

    // must be deleted by the caller method
    vector<long> chk = vector<long>(num_vertices, -1);

    long component = 0;
    for (long u=0; u<num_vertices; ++u)
    {
        if (chk[u] < 0)
            dfs_to_tag_components(u, component, chk, adj_list);

        ++component;
    }

    return chk;
}

///////////////////////////////////////////////////////////////////////////////

KStabCutGenerator::KStabCutGenerator(GRBModel *model, GRBVar *x_vars, IO *instance)
{
    this->model = model;
    this->x_vars = x_vars;
    this->instance = instance;
    this->num_vars = instance->graph->num_edges;

    this->oci_counter = 0;
    this->oci_stats = new oci_statistics();
}

KStabCutGenerator::~KStabCutGenerator()
{
    delete oci_stats;
}

void KStabCutGenerator::callback()
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

            // generate cuts only at root node?
            if (CUTS_AT_ROOT_ONLY && getDoubleInfo(GRB_CB_MIPNODE_NODCNT) > 0)
                return;

            if (OCI_STRATEGY_SHIFT_AFTER_ROOT)
                if (getDoubleInfo(GRB_CB_MIPNODE_NODCNT) > 0)
                    OCI_STRATEGY = OCI_STRATEGY_SHIFT_CHOICE;

            x_val = this->getNodeRel(x_vars, num_vars);

            // find violated odd-cycle inequalities (if any) and add cuts to the model
            if (SEPARATE_OCI)
                run_oci_separation(ADD_USER_CUTS);

            delete[] x_val;
        }

        // callback from a new MIP incumbent: including LAZY CONSTRAINTS (NOT THE CASE FOR OCI)
        /*
        else if (where == GRB_CB_MIPSOL)
        {
            x_val = this->getSolution(x_vars, num_vars);

            // find violated odd-cycle inequalities (if any) and add cuts to the model
            if (SEPARATE_OCI)
                run_oci_separation(ADD_LAZY_CNTRS);

            delete[] x_val;
        }
        */
    }
    catch (GRBException e)
    {
        cout << "Error " << e.getErrorCode() << " during KStabCutGenerator::callback(): ";
        cout << e.getMessage() << endl;
    }
    catch (...)
    {
        cout << "Unexpected error during KStabCutGenerator::callback()" << endl;
    }
}


bool KStabCutGenerator::separate_lpr()
{
    /// Interface to be used when solving the LP relaxation only.

    try
    {
        bool model_updated = false;
        
        x_val = new double[num_vars];
        for (long i=0; i < num_vars; ++i)
            x_val[i] = x_vars[i].get(GRB_DoubleAttr_X);

        // find violated oci's (if any) and add cuts to the model
        if (SEPARATE_OCI)
            model_updated = run_oci_separation(ADD_STD_CNTRS);

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

bool KStabCutGenerator::run_oci_separation(int kind_of_cuts)
{
    /// wrapper for the separation procedure to suit different kinds of cuts

    bool model_updated = false;

    // eventual cuts are stored here
    vector<GRBLinExpr> cuts_lhs = vector<GRBLinExpr>();
    vector<long> cuts_rhs = vector<long>();

    // run separation algorithm by Gerards & Schrijver (1986)
    model_updated = separate_oci(cuts_lhs,cuts_rhs);

    if (model_updated)
    {
        // add cuts
        for (unsigned long idx = 0; idx<cuts_lhs.size(); ++idx)
        {
            ++oci_counter;

            if (kind_of_cuts == ADD_USER_CUTS)
                addCut(cuts_lhs[idx] <= cuts_rhs[idx]);

            /*
            // OCIs are not used as lazy constraints - only user cuts
            else if (kind_of_cuts == ADD_LAZY_CNTRS)
                addLazy(cuts_lhs[idx] <= cuts_rhs[idx]);
            */

            else // kind_of_cuts == ADD_STD_CNTRS
                model->addConstr(cuts_lhs[idx] <= cuts_rhs[idx]);
        }
    }

    return model_updated;
}

bool KStabCutGenerator::separate_oci( vector<GRBLinExpr> &cuts_lhs,
                                      vector<long> &cuts_rhs )
{
    /***
     * Solve the separation problem for odd-cycle inequalities, assuming that
     * the current point satisfies the trivial edge inequalities.
     * Original reference: Gerards & Schrijver (1986)
     */

    // prevent floating point errors by ignoring digits beyond set precision 
    for (long i=0; i < this->num_vars; ++i)
    {
        double tmp = x_val[i] * std::pow(10,OCI_SEPARATION_PRECISION_IN_IP);
        tmp = std::round(tmp);
        x_val[i] = tmp * std::pow(10,-OCI_SEPARATION_PRECISION_IN_IP);
    }

    long separated = 0;

    // 1. check vars equal to 0 or 1
    long reduce_count = 0;
    vector<bool> reduce = vector<bool>(num_vars, false);

    for (long i=0; i<num_vars; ++i)
    {
        if (x_val[i] <= EPSILON_TOL || x_val[i] >= 1-EPSILON_TOL)
        {
            reduce.at(i) = true;
            ++reduce_count;
        }
    }

    // 2. calculate weights for C (edges in the conflict graph of an MSTCC instance)

    vector<double> c_weight;

    /***
     * Elementary sparse matrix structure to avoid storing an explicit adjacency
     * matrix (prohibitive for large instances) to get weight of an edge from
     * its endpoints' indices.
     * Also, to avoid calling the expensive insert operation on the map
     * redundantly, we only store one element for each adjacency: the map to 
     * use is that of the smaller index, and the map key is the larger one.
     */
    vector< map<long,double> > c_weight_from_vertices;
    for (long i=0; i<num_vars; ++i)
        c_weight_from_vertices.push_back(map<long,double>());

    for (long p=0; p<instance->num_conflicts; ++p)
    {
        pair<long,long> edges = instance->conflicts.at(p);
        long e1 = edges.first;
        long e2 = edges.second;
        
        // weight c(i,j) = (1 - y_i - y_j)/2
        double w = 1.0;
        w -= x_val[e1];
        w -= x_val[e2];
        w /= 2.0;
        
        // avoid null-cost edges (implementation tweak #2 in rebennack et al., 2012)
        c_weight.push_back(w <= EPSILON_TOL ? EPSILON_TOL : w);
        
        if (e1 < e2)
            c_weight_from_vertices[e1][e2] = c_weight.back();
        else
            c_weight_from_vertices[e2][e1] = c_weight.back();
    }

    // 3. adjacency list of the auxiliary bipartite graph H (doubling vertices and edges)
    vector< vector<long> > aux_adj_list;
    long aux_vertex_count = 0;
    long aux_arc_count = 0;

    // handle to map a vertex of H to the original vertex of G
    vector<long> map_h2g = vector<long>();

    // handle to map an original vertex 'e' of G to the first vertex 'e+' of H
    vector<long> map_g2h = vector<long>();

    // vertices
    for (long i=0; i<num_vars; ++i)
    {
        // tweak #1 in rebennack et al. (2012): no need to consider in H vars equal to 0 or 1
        if (reduce.at(i))
            map_g2h.push_back(-1);  // e is not included in H
        else
        {
            // two vertices (e+ and e-) in H for each vertex e in G
            aux_adj_list.push_back(vector<long>());
            aux_adj_list.push_back(vector<long>());

            map_g2h.push_back(aux_vertex_count);

            map_h2g.push_back(i);
            map_h2g.push_back(i);

            aux_vertex_count += 2;
        }
    }

    // arcs
    for (vector<pair<long,long> >::iterator it = instance->conflicts.begin(); it != instance->conflicts.end(); ++it)
    {
        long v1 = (*it).first;
        long v2 = (*it).second;

        if (!reduce.at(v1) && !reduce.at(v2))
        {
            long v1a = map_g2h[v1];
            long v1b = v1a + 1;

            long v2a = map_g2h[v2];
            long v2b = v2a + 1;

            // arcs between v1+ and v2-
            aux_adj_list[v1a].push_back(v2b);
            aux_adj_list[v2b].push_back(v1a);

            // arcs between v1- and v2+
            aux_adj_list[v1b].push_back(v2a);
            aux_adj_list[v2a].push_back(v1b);

            aux_arc_count += 4;
        }
    }

    // auxiliary structures for storing and selecting cuts to add
    double most_violated = 0;
    long most_violated_idx = -1;
    vector<violated_oci*> cuts;

    // get connected components to avoid computing sp for unconnected vertices
    vector<long> components = check_components(aux_adj_list);

    // setup complete, starting search from each source vertex
    for (long source=0; source<num_vars; ++source)
    {
        if (!reduce.at(source))
        {
            long u1 = map_g2h[source];
            long u2 = u1+1;

            // TODO: REMOVE ME
            if (u1 < 0)
                cerr << endl << "SHOULD NEVER OCCUR: tried to index a non-reduced vertex before solving sssp" << endl << endl;

            // source and destination are not in different components
            if (components.at(u1) == components.at(u2))
            {
                // 4. single-source shortest-path from u1 to u2 using Dijkstra's algorithm

                // initialize vertices (heap elements) and priority queue
                vector<sssp_heap_element*> heap_nodes = vector<sssp_heap_element*>();

                for (long i=0; i<aux_vertex_count; ++i)
                {
                    sssp_heap_element *handle = new sssp_heap_element();
                    handle->estimate = numeric_limits<double>::max();
                    handle->v = i;
                    handle->previous = -1;
                    handle->heap_pos = i+1; // head implementation has a dummy vertex

                    heap_nodes.push_back(handle);
                }

                heap_nodes[u1]->estimate = 0;  // source

                sssp_minheap queue;
                queue.build_min_heap(heap_nodes);

                // kernel: iteratively select closest vertex, and relax incident edges
                bool found_u2 = false;
                while (queue.get_size() > 0 && !found_u2)
                {
                    sssp_heap_element *i = queue.extract_min();

                    // NOTE: in our case, we are able to stop as soon as goal vertex 'u2' is chosen
                    if (i->v == u2)
                        found_u2 = true;
                    else
                        for (long adj=0; (unsigned) adj < aux_adj_list[i->v].size(); ++adj)
                        {
                            // arc (i,j) information
                            sssp_heap_element *j = heap_nodes[aux_adj_list[i->v][adj]];

                            long i_original = map_h2g[ i->v ];
                            long j_original = map_h2g[ j->v ];

                            double w = (i_original < j_original) ?
                                c_weight_from_vertices[i_original][j_original] :
                                c_weight_from_vertices[j_original][i_original] ;

                            // relax arc(i,j)
                            if (j->estimate > (double) i->estimate + w)
                            {
                                j->estimate = (double) i->estimate + w;
                                j->previous = i->v;

                                // update heap
                                if (queue.decrease_key(j->heap_pos, i->estimate + w) == false)
                                    cout << "UNEXPECTED ERROR IN PRIORITY QUEUE (DECREASING KEY WITH HIGHER VALUE)" << endl;
                            }
                        }
                }

                // retrieve shortest path
                vector<long> shortest_path = vector<long>();
                shortest_path.push_back(map_h2g[u2]);   // corresponding vertex in G
                long aux = u2;
                while (heap_nodes[aux]->previous >= 0) //!= u1)
                {
                    long prev = heap_nodes[aux]->previous;
                    shortest_path.push_back(map_h2g[prev]);
                    aux = prev;
                }

                // 5. get odd cycle in G from walk corresponding to shortest path in H
                vector<bool> chk_v = vector<bool>(num_vars, false);

                chk_v[shortest_path[0]] = true;   // first vertex start labeled as visited

                /***
                 * Just like with c_weight_from_vertices, we have to avoid an
                 * explicit quadractic structure for checking/labeling edges.
                 */
                vector< map<long,bool> > chk_e;
                for (long i=0; i<num_vars; ++i)
                    chk_e.push_back(map<long,bool>());

                vector<vector<long> > cycles_found = vector<vector<long> >();

                bool done = false;
                unsigned long pointer = 0;
                while (!done)
                {
                    long u = shortest_path[pointer];
                    long v = shortest_path[pointer+1];

                    // check for repeated vertex
                    if (chk_v[v] == true)
                    {
                        // found a cycle
                        vector<long> cycle = vector<long>();

                        // store vertices in this cycle and remove from the path
                        cycle.push_back(v);
                        cycle.push_back(u);

                        shortest_path.erase(shortest_path.begin()+pointer+1);  // iterator to v
                        shortest_path.erase(shortest_path.begin()+pointer);    // iterator to u
                        chk_v[u] = 0;

                        if (u <= v)
                            chk_e[u][v] = false;
                        else
                            chk_e[v][u] = false;

                        long terminus = v;   // complete cycle
                        --pointer;
                        while (shortest_path[pointer] != terminus)
                        {
                            v = u;
                            u = shortest_path[pointer];
                            cycle.push_back(u);

                            shortest_path.erase(shortest_path.begin()+pointer);  // iterator to pointer
                            chk_v[u] = 0;

                            if (u <= v)
                                chk_e[u][v] = false;
                            else
                                chk_e[v][u] = false;

                            --pointer;
                        }

                        cycles_found.push_back(cycle);

                        // check if we're done
                        if (pointer == shortest_path.size() - 1)
                            done = true;
                        else
                        {
                            u = shortest_path[pointer];
                            v = shortest_path[pointer+1];

                            long smaller_idx = (u<=v) ? u : v;
                            long greater_idx = (smaller_idx==u) ? v : u;

                            // check for repeated edges
                            while ( !done && 
                                    chk_e[smaller_idx].count(greater_idx) > 0 &&  // element is there
                                    chk_e[smaller_idx][greater_idx] )             // and it is 1
                            {
                                shortest_path.erase(shortest_path.begin()+pointer+1);  // iterator to v
                                shortest_path.erase(shortest_path.begin()+pointer);    // iterator to u
                                chk_v[u] = 0;

                                if (u <= v)
                                    chk_e[u][v] = false;
                                else
                                    chk_e[v][u] = false;

                                --pointer;
                                u = shortest_path[pointer];

                                // check if we're done
                                if (pointer == shortest_path.size() - 1)
                                    done = true;
                                else
                                    v = shortest_path[pointer+1];
                            }
                        }
                    }
                    else  // new arc: label new vertex and arcs as visited and proceed
                    {
                        chk_v[v] = true;

                        if (u <= v)
                            chk_e[u][v] = true;
                        else
                            chk_e[v][u] = true;

                        ++pointer;
                    }
                }

                long num_cycles = cycles_found.size();

                // 6. CHECK EACH ODD CYCLE FOR VIOLATION OF THE CORRESPONDING CONSTRAINT
                for (long i=0; i < num_cycles; ++i)
                {
                    vector<long> &cycle = cycles_found[i];
                    long len = cycle.size();

                    // check odd cycles only
                    if (len % 2 == 1)
                    {
                        double cost = 0.;

                        for (long j=0; j<len; ++j)
                        {
                            long v1 = cycle[j];
                            long v2 = cycle[(j+1)%len];

                            if (v1 < v2)
                                cost += (double) c_weight_from_vertices[v1][v2];
                            else
                                cost += (double) c_weight_from_vertices[v2][v1];
                        }

                        // check if corresponding cut is violated
                        if (cost < 0.5 - OCI_VIOLATION_TOL_IN_IP)
                        {
                            ++separated;

                            // cut object
                            violated_oci *oci = new violated_oci(num_vars);
                            double current_sum = 0.;

                            for (long j=0; j<len; ++j)
                            {
                                long v = cycle[j];
                                oci->C.push_back(v);
                                oci->coefficients[v] = 1;

                                current_sum += x_val[v];
                            }

                            double tmp = (double) (len - 1.)/2.;
                            oci->infeasibility = (double) current_sum - tmp;

                            // 7. ADD CUT
                            if (OCI_STRATEGY == ALL_CUTS)
                            {
                                oci_stats->oci_len[len] = oci_stats->oci_len[len] <= 0 ?
                                                          1 : oci_stats->oci_len[len] + 1;

                                if (STORE_OCI_CUT_POOL)
                                    oci_stats->oci_pool[oci->toString()] = 1;

                                // store oci (caller method adds it to the model)
                                GRBLinExpr violated_constraint = 0;
                                for (long j=0; j<len; ++j)
                                    violated_constraint += x_vars[cycle[j]];

                                cuts_lhs.push_back(violated_constraint);
                                cuts_rhs.push_back((double) (len-1.) / 2.);

                                delete oci;
                            }
                            else
                            {
                                // only store oci object, will select cuts at the end
                                cuts.push_back(oci);
 
                                if (oci->infeasibility > most_violated)
                                {
                                    most_violated = oci->infeasibility;
                                    most_violated_idx = cuts.size() - 1;
                                }
                            }

                        } // violated cut test
                    }

                } // done checking this cycle

                // clean up
                for (vector<sssp_heap_element*>::iterator it = heap_nodes.begin(); it != heap_nodes.end(); ++it)
                    delete *it;
                heap_nodes.clear();
            }
        }

    } // done with this source vertex

    // 8. ENHANCED CUT ADDITION STRATEGIES

    if (separated > 0 && OCI_STRATEGY != ALL_CUTS)
    {
        // 8.1 ADD MOST VIOLATED CUT, INCLUDED IN BOTH ENHANCED STRATEGIES

        long len = cuts[most_violated_idx]->C.size();

        oci_stats->oci_len[len] = oci_stats->oci_len[len] <= 0 ?
                                  1 : oci_stats->oci_len[len] + 1;

        if (STORE_OCI_CUT_POOL)
            oci_stats->oci_pool[cuts[most_violated_idx]->toString()] = 1;

        // lhs
        GRBLinExpr violated_constraint = 0;
        for (long j=0; j<len; ++j)
        {
            long node_in_this_cycle = cuts[most_violated_idx]->C.at(j);
            violated_constraint += x_vars[node_in_this_cycle];
        }

        // store oci (caller method adds it to the model)
        cuts_lhs.push_back(violated_constraint);
        cuts_rhs.push_back((double) (len-1.) / 2.);

        // we are done if the strategy is to add just the most violated cut

        if (OCI_STRATEGY == ORTHOGONAL_CUTS)
        {
            /***
             * 6.2 ADD ANY OTHER CUT WHOSE HYPERPLANE HAS INNER PRODUCT WITH 
             * THAT OF THE MOST VIOLATED CUT CLOSE TO ZERO
             */

            // 2-norm of the vector corresponding to most violated inequality
            double norm_v1 = 0.;
            for(long i=0; i<num_vars; ++i)
            {
                double base = (double) cuts[most_violated_idx]->coefficients[i];
                double sq = pow(base, 2.);
                norm_v1 += sq;
            }
            //double rhs1 = ((double) (cuts[most_violated_idx]->C.size()-1.) / 2.);
            //norm_v1 += pow(rhs1, 2);
            norm_v1 = sqrt(norm_v1);

            for (long idx=0; (unsigned) idx<cuts.size(); ++idx)
            {
                if (idx != most_violated_idx)
                {
                    // norm-2 of candidate cut vector
                    double norm_v2 = 0.;
                    for(long i=0; i<num_vars; ++i)
                    {
                        double base = (double) cuts[idx]->coefficients[i];
                        double sq = pow(base, 2.);
                        norm_v2 += sq;
                    }
                    //double rhs2 = ((double) (cuts[idx]->C.size()-1.) / 2.);
                    //norm_v2 += pow(rhs2, 2);
                    norm_v2 = sqrt(norm_v2);

                    // inner product
                    double dot = 0.;
                    for(long i=0; i<num_vars; ++i)
                    {
                        double v1_i = (double) cuts[most_violated_idx]->coefficients[i];
                        double v2_i = (double) cuts[idx]->coefficients[i];
                        dot += v1_i * v2_i;
                    }
                    //dot += (rhs1*rhs2);  // adds because both rhs are non-negative

                    // add cut if normalized dot product is close to 0
                    double norm = norm_v1 * norm_v2;
                    double inner_prod = (double) dot / norm;

                    if (inner_prod <= OCI_ORTHOGONALITY_TOL)
                    {
                        long len = cuts[idx]->C.size();

                        oci_stats->oci_len[len] = oci_stats->oci_len[len] <= 0 ?
                                                  1 : oci_stats->oci_len[len] + 1;

                        if (STORE_OCI_CUT_POOL)
                            oci_stats->oci_pool[cuts[idx]->toString()] = 1;

                        // lhs
                        GRBLinExpr violated_constraint = 0;
                        for (long j=0; j<len; ++j)
                        {
                            long node_in_this_cycle = cuts[idx]->C.at(j);
                            violated_constraint += x_vars[node_in_this_cycle];
                        }

                        // store oci (caller method adds it to the model)
                        cuts_lhs.push_back(violated_constraint);
                        cuts_rhs.push_back((double) (len-1.) / 2.);
                    }
                }

            } // checked each cut
        }

    } // done with cut addition in enhanced strategies

    // clean up
    for (vector<violated_oci*>::iterator it = cuts.begin(); it != cuts.end(); ++it)
        delete *it;
    cuts.clear();

    return (separated > 0);
}
