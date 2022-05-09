#ifndef _KSTAB_CUT_GEN_H_
#define _KSTAB_CUT_GEN_H_

#include <iostream>
#include <sstream>
#include <iomanip>

#include "gurobi_c++.h"

#include "io.h"
#include "kstab_model.h"

// cut selection strategies
#define MOST_VIOLATED_CUT 1
#define ORTHOGONAL_CUTS 2
#define ALL_CUTS 3

// kinds of cuts
#define ADD_USER_CUTS 1
#define ADD_LAZY_CNTRS 2
#define ADD_STD_CNTRS 3

/***
 * \file kstab_cut_generator.h
 * 
 * Module for the Gurobi callback class, implementing the separation procedure
 * for odd-cycle inequalities (OCI).
 * 
 * Extends (and implements) the abstract base class in Gurobi.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 06.05.2022
 */

class oci_statistics;

class KStabCutGenerator: public GRBCallback
{
public:
    KStabCutGenerator(GRBModel *, GRBVar*, IO*);
    virtual ~KStabCutGenerator();

protected:
    friend class KStabModel;
    friend class StableSpanningTreeModel;
    friend class violated_oci;

    void callback();
    bool separate_lpr();

    // input and model data
    IO *instance;
    GRBModel *model;

    GRBVar* x_vars;
    double *x_val;
    long num_vars;

    long oci_counter;
    oci_statistics *oci_stats;
    bool run_oci_separation(int);
    bool separate_oci(vector<GRBLinExpr> &, vector<long> &);
};

class oci_statistics
{
public:
    map<long,long> oci_len;             // length of odd cycles found
    map<string,long> oci_pool;          // optionally store cuts
};

/// information of a violated subtour elimination constraint
class violated_oci
{
public:
    violated_oci(long vertex_count)
    {
        this->vertex_count = vertex_count;
        coefficients = vector<bool>(vertex_count, false);
    }

    virtual ~violated_oci() { }
    
    string toString()
    {
        stringstream oci_lhs;
        oci_lhs.str("");
        for (long i = 0; i<vertex_count; ++i)
            oci_lhs << coefficients.at(i);
        return oci_lhs.str();
    }

    long vertex_count;          // instance parameter
    vector<long> C;             // index of vertices in a cycle
    double infeasibility;       // amount by which the oci is violated
    vector<bool> coefficients;  // hyperplane coefficients
};

#endif
