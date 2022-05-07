#include "kstab_cut_generator.h"

/// algorithm setup switches

bool CUTS_AT_ROOT_ONLY = false;

int  OCI_STRATEGY = ALL_CUTS;
bool STORE_OCI_CUT_POOL = false;

#define OCI_ORTHOGONALITY_TOL 0.1
#define OCI_VIOLATION_TOL_IN_IP 0.0001
#define OCI_SEPARATION_PRECISION_IN_IP 14   // <= 14 without changing everything to long double (which gurobi ignores)

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
            if (CUTS_AT_ROOT_ONLY && getDoubleInfo(GRB_CB_MIPNODE_NODCNT) != 0.)
                return;

            x_val = this->getNodeRel(x_vars, num_vars);

            // find violated odd-cycle inequalities (if any) and add cuts to the model
            run_oci_separation(ADD_USER_CUTS);

            delete[] x_val;
        }

        // callback from a new MIP incumbent: including LAZY CONSTRAINTS (NOT THE CASE FOR OCI)
        /*
        else if (where == GRB_CB_MIPSOL)
        {
            x_val = this->getSolution(x_vars, num_vars);

            // find violated odd-cycle inequalities (if any) and add cuts to the model
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
        x_val = new double[num_vars];
        for (long i=0; i < num_vars; ++i)
            x_val[i] = x_vars[i].get(GRB_DoubleAttr_X);

        // find violated subtour eliminaton constraints (if any) and add cuts to the model
        bool model_updated = run_oci_separation(ADD_STD_CNTRS);

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
     * Solve the separation problem for odd-cycle inequalities.
     * Original reference: Gerards & Schrijver (1986)
     */

    // prevent floating point errors by ignoring digits beyond set precision 
    for (long i=0; i < this->num_vars; ++i)
    {
        double tmp = x_val[i] * std::pow(10,OCI_SEPARATION_PRECISION_IN_IP);
        tmp = std::round(tmp);
        x_val[i] = tmp * std::pow(10,-OCI_SEPARATION_PRECISION_IN_IP);
    }

    // TO DO: ALL

    return false;
}
