#include "mstcc_model.h"

// cut selection strategies
#define FIRST_CUT 0
#define MOST_VIOLATED_CUT 1
#define ORTHOGONAL_CUTS 2
#define ALL_CUTS 3

#define ORTHOGONALITY_TOL 0.1

// algorithm setup switches
int USE_SEC_STRATEGY = ORTHOGONAL_CUTS;
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
}


StableSpanningTreeModel::~StableSpanningTreeModel()
{
    // TODO: nothing here?
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
            
            for (long i=0; i < this->instance->graph->num_edges; ++i)
                cout << "x[" << i << "] = " << this->x[i].get(GRB_DoubleAttr_X) << endl;
            #endif

            // eventual cuts are stored here
            vector<GRBLinExpr> cuts_lhs;
            vector<long> cuts_rhs;

            if (USE_FAST_INTEGER_CUT)
            {
                // more efficient separation procedure for an integral solution
                model_updated = separate_SEC_integer(cuts_lhs,cuts_rhs);
            }

            // if not using the integer separation or if solution is fractional
            if (!model_updated)
            {
                // standard separation procedure
                model_updated = separate_SEC(cuts_lhs,cuts_rhs);
            }

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


void StableSpanningTreeModel::dfs_to_count( long u, 
                                            bool *chk_v, 
                                            bool **chk_e, 
                                            vector<long> &seen_edges, 
                                            vector< vector<long> > &adj_list, 
                                            long *v_count, 
                                            long *e_count )
{
    /// specialized depth-first search to count vertices and edges in subgraph

    chk_v[u] = 1;
    ++(*v_count);
    
    for (unsigned i=0; i<adj_list[u].size(); ++i)
    {
        int v = adj_list[u].at(i);
        
        if (chk_e[u][v] == 0 && chk_e[v][u] == 0)
        {
            // every edge we see is counted (e.g. whether back edge or not)
            ++(*e_count);

            // sets only (u,v) to find violated constraint more easily
            chk_e[u][v] = 1;
            seen_edges.push_back( this->instance->graph->index_matrix[u][v] );
        }
        
        if (chk_v[v] == 0)
            dfs_to_count(v, chk_v, chk_e, seen_edges, adj_list, v_count, e_count);
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

    #ifdef DEBUG
    cout << "trying fast integer cut... ";
    #endif

    // proceed only if the current solution is integer
    if (s_.size() == (unsigned) instance->graph->num_edges)
    {
        // makes adjlist of the subgraph induced by the current solution
        vector< vector<long> > s_adj_list;
        for (long i=0; i<instance->graph->num_vertices; ++i)
            s_adj_list.push_back(vector<long>());
        
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
        
        // check if solution contains a cycle
        bool *chk_v = new bool[instance->graph->num_vertices];
        memset(chk_v, 0, sizeof(bool)*instance->graph->num_vertices);
        
        for (long root=0; root<instance->graph->num_vertices; ++root)
        {
            if (!chk_v[root])
            {
                // checks vertices and edges as visited
                long count_v = 0;
                long count_e = 0;
                vector<long> list_e; vector<long>();
                bool **chk_e = new bool*[instance->graph->num_vertices];
                for (long i=0; i<instance->graph->num_vertices; ++i)
                {
                    chk_e[i] = new bool[instance->graph->num_vertices];
                    memset(chk_e[i], 0, sizeof(bool)*instance->graph->num_vertices);
                }
                
                dfs_to_count(root, chk_v, chk_e, list_e, s_adj_list, &count_v, &count_e);
                
                // a component includes a cycle if #edges >= vertices
                if (count_e >= count_v)
                {
                    // set S with checked vertices yields violated constraint
                    GRBLinExpr violated_constraint = 0;
                    
                    #ifdef DEBUG
                    cout << "succeeded on ";
                    #endif

                    // lhs
                    for (vector<long>::iterator it = list_e.begin(); it != list_e.end(); ++it)
                    {
                        violated_constraint += x[*it];

                        #ifdef DEBUG
                        cout << "+ x[" << *it << "]";
                        #endif
                    }

                    // save this SEC in the reference arg (added selectively to model by the caller function)
                    cuts_lhs.push_back(violated_constraint);
                    cuts_rhs.push_back(count_v-1);

                    #ifdef DEBUG
                    cout << " <= " << count_v-1 << endl;
                    #endif

                    // clean up
                    for (long i = 0; i<instance->graph->num_vertices; ++i)
                        delete[] chk_e[i];
                    delete[] chk_e;
                    delete[] chk_v;

                    return true;
                }
                
                for (long i = 0; i<instance->graph->num_vertices; ++i)
                    delete[] chk_e[i];
                delete[] chk_e;
            }
        }
        
        delete[] chk_v;
    }

    #ifdef DEBUG
    cout << "failed" << endl;
    #endif
    return false;
}


bool StableSpanningTreeModel::separate_SEC( vector<GRBLinExpr> &cuts_lhs,
                                            vector<long> &cuts_rhs )
{
    /***
     * solve the separation problem for subtour elimination constraints
     */

    cout << endl << "will try to separate SEC here..." << endl << endl;
    return cuts_lhs.size() != cuts_rhs.size(); //==false
}