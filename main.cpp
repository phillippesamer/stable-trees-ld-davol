/***
 * \file main.cpp
 * 
 * Lagrangean Decomposition based dual ascent algorithm to compute lower bounds
 * for minimum weight stable (or conflict-free) spanning trees.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 01.11.2021
 */

#include "io.h"
#include "kstab_model.h"
#include "mstcc_model.h"
#include "ldda.h"

#include <cstdlib>
#include <fstream>

using namespace std;

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cout << "usage: \t" << argv[0] << " [input file]" << endl << endl;
        return 0;
    }

    IO* instance = new IO();
    
    // parse given input file and look for errors in it
    if (instance->parse_gcclib(argv[1]) == false)
    {
        cout << "unable to parse input file" << endl;
        delete instance;
        return 0;
    }

    // LP relaxation bound of the MSTCC natural IP formulation
    StableSpanningTreeModel *lpr_model = new StableSpanningTreeModel(instance);

    if ( lpr_model->solve_lp_relax(false) )
    {
        KStabModel *model = new KStabModel(instance);
        LDDA *lagrangean = new LDDA(instance, model);

        // trivial combinatorial bounds from min weight spanning tree and kstab
        instance->run_mst();
        model->solve(true);

        // TO DO: remove this (temporary output for experiments)
        stringstream table_row;

        cout << "_____________________________________________________________________________" << endl << endl;

        cout << "kstab bound: " << model->solution_weight
             << " (runtime " << fixed << model->solution_runtime << ")" << endl;

        cout << "mst bound: " << instance->get_mst_weight()
             << " (runtime " << fixed << instance->get_mst_runtime() << ")" << endl;

        cout << "lp bound: " << lpr_model->lp_bound
             << " (" << lpr_model->lp_passes << " passes,"
             << " runtime: " << fixed << lpr_model->lp_runtime << ")" << endl;
        cout << "_____________________________________________________________________________" << endl << endl;

        // TO DO: remove this (temporary output for experiments)
        table_row << left;
        table_row << setw(50) << argv[1];
        table_row << setw(5) << "  &  ";
        table_row << setw(10) << model->solution_weight;
        table_row << setw(5) << "  &  ";
        table_row << setw(10) << fixed << model->solution_runtime;
        table_row << setw(5) << "  &  ";
        table_row << setw(10) << instance->get_mst_weight();
        table_row << setw(5) << "  &  ";
        table_row << setw(10) << fixed << instance->get_mst_runtime();
        table_row << setw(5) << "  &  ";
        table_row << setw(10) << lpr_model->lp_bound;
        table_row << setw(5) << "  &  ";
        table_row << setw(10) << fixed << lpr_model->lp_runtime;
        table_row << setw(5) << "  &  ";

        delete lpr_model;

        // Lagrangean Decomposition bound
        if ( lagrangean->dual_ascent(false) == false)
        {
            // TO DO: remove this (temporary for experiments)
            lagrangean->bound_log.push_back(numeric_limits<long>::min());
        }

        cout << endl
             << "ldda bound: " << lagrangean->bound_log.back()
             << " (runtime " << fixed << lagrangean->runtime << ")" << endl;

        // write log file (input file name + "_ldda.log")
        stringstream log = lagrangean->create_log();

        char buffer[200];
        int cx = snprintf(buffer, 200, "%s_ldda.log", argv[1]);
        ofstream logfile(buffer);
        if (cx>=0 && cx<200 && logfile.is_open())
        {
            logfile << log.str();
            logfile.close();
        }
        else
        {
            cout << log.str();
            cout << "ERROR: unable to write log file; dumped to screen" << endl;
        }

        // TO DO: remove this (temporary output for experiments)
        table_row << setw(10) << lagrangean->bound_log.back();
        table_row << setw(5) << "  &  ";
        table_row << setw(10) << fixed << lagrangean->runtime;
        table_row << setw(5) << endl;

        ofstream xpfile("xp1table.dat", ofstream::app);
        if (xpfile.is_open())
        {
            xpfile << table_row.str();
            xpfile.close();
        }
        else
        {
            cout << "ERROR: unable to write dat file; dumping to screen:" << endl;
            cout << table_row.str();
        }

        // clean up
        delete lagrangean;
        delete model;
    }
    else
    {
        delete lpr_model;
    }

    delete instance;
    return 0;
}
