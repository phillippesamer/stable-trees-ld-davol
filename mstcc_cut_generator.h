#ifndef _MSTCC_CUT_GEN_H_
#define _MSTCC_CUT_GEN_H_

#include <iostream>
#include <sstream>
#include <iomanip>

#include "gurobi_c++.h"

#include "io.h"
#include "kstab_model.h"
#include "kstab_cut_generator.h"
#include "mstcc_model.h"

/***
 * \file mstcc_cut_generator.h
 * 
 * Module for the Gurobi callback class, implementing the separation procedure
 * for subtour elimination constraints (SEC).
 * 
 * Extends KStabCutGenerator, the corresponding cut generating class for 
 * KStabModel (base class of StableSpanningTreeModel).
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 29.03.2022
 */

class sec_statistics;

class SSTCutGenerator: public KStabCutGenerator
{
public:
    SSTCutGenerator(GRBModel *, GRBVar*, IO*);
    virtual ~SSTCutGenerator();

protected:
    friend class KStabModel;
    friend class StableSpanningTreeModel;
    friend class violated_sec;

    void callback();
    bool separate_lpr();

    long sec_counter;
    sec_statistics *sec_stats;
    bool run_sec_separation(int);
    bool separate_sec_integer(vector<GRBLinExpr> &, vector<long> &);
    bool separate_sec_folklore(vector<GRBLinExpr> &, vector<long> &);
    bool separate_sec_classical(vector<GRBLinExpr> &, vector<long> &);
    void dfs_checking_acyclic(long, long, vector<bool> &, long &, stack<long> &, vector< list<long> > &, bool &);
};

class sec_statistics
{
public:
    map<string,long> pool;               // optionally store cuts

    vector<long> sec_diff_cuts;          // # of different cuts in each separation
    vector<double> sec_infeas_std_dev;   // std.dev. of the amount violated    
};

/// information of a violated subtour elimination constraint
class violated_sec
{
public:
    violated_sec(long edge_count)
    {
        this->edge_count = edge_count;
        coefficients = vector<bool>(edge_count, false);
    }

    virtual ~violated_sec() { }
    
    string toString()
    {
        stringstream sec_lhs;
        sec_lhs.str("");
        for (long i=0; i<edge_count; ++i)
            sec_lhs << coefficients.at(i);
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

bool inline sort_pairs_by_snd_val(pair<long,double> a, pair<long,double> b)
{
    return ( a.second < b.second ); 
}

#endif
