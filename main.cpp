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

// execution switches
double RUN_KSTAB_WITH_TIME_LIMIT = 3600;

bool RUN_STEEPEST_ASCENT_LDDA = false;
bool WRITE_LDDA_LOG_FILE = true;

bool APPEND_SUMMARY_TO_DAT_FILE = true;
string SUMMARY_FILE_NAME = string("xp1table.dat");

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
        // trivial combinatorial bounds from min weight spanning tree and kstab
        instance->run_mst();

        KStabModel *model = new KStabModel(instance);
        model->set_time_limit(RUN_KSTAB_WITH_TIME_LIMIT);
        model->solve(true);

        cout << "_____________________________________________________________________________" << endl << endl;

        cout << "kstab bound: " << model->solution_dualbound
             << " (runtime " << fixed << model->solution_runtime << ")";
        if (model->solution_status != AT_OPTIMUM)
            cout << " *** NOT OPTIMAL ***";
        cout << endl;

        cout << "mst bound: " << instance->get_mst_weight()
             << " (runtime " << fixed << instance->get_mst_runtime() << ")" << endl;

        cout << "lp bound: " << lpr_model->lp_bound
             << " (" << lpr_model->lp_passes << " passes,"
             << " runtime " << fixed << lpr_model->lp_runtime << ")" << endl;
        cout << "_____________________________________________________________________________" << endl << endl;

        // summary output file (for experiments with many instances)
        stringstream table_row;
        if (APPEND_SUMMARY_TO_DAT_FILE)
        {
            table_row << left;
            table_row << setw(50) << argv[1];
            table_row << setw(5) << "  &  ";
            if (model->solution_status == AT_OPTIMUM)
                table_row << setw(25) << model->solution_weight;
            else if (model->solution_status == IS_INFEASIBLE)
                table_row << setw(25) << "x";
            else
            {
                stringstream tmp_str;
                tmp_str << model->solution_weight << " ?";
                table_row << setw(25) << tmp_str.str();
            }
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
        }

        delete lpr_model;

        if (model->solution_status == AT_OPTIMUM)
        {
            // Lagrangean Decomposition bound
            LDDA *lagrangean = new LDDA(instance, model);
            bool ldda_complete = lagrangean->dual_ascent(RUN_STEEPEST_ASCENT_LDDA);

            cout << endl << "ldda bound: ";
            if ( ldda_complete)
                cout << lagrangean->bound_log.back();
            else
                cout << " - ";
            cout << " (runtime " << fixed << lagrangean->runtime << ")" << endl;

            if (WRITE_LDDA_LOG_FILE)
            {
                // write LDDA log file (input file name + "_ldda.log")
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
            }

            if (APPEND_SUMMARY_TO_DAT_FILE)
            {
                if (ldda_complete)
                    table_row << setw(10) << lagrangean->bound_log.back();
                else
                    table_row << setw(10) << " - ";
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << fixed << lagrangean->runtime;
                table_row << setw(5) << endl;
            }

            delete lagrangean;
        }
        else
        {
            // kstab model could not be solved to optimality within time limit

            if (APPEND_SUMMARY_TO_DAT_FILE)
            {
                table_row << setw(10) << " - ";
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << " - ";
                table_row << setw(5) << endl;
            }
        }

        if (APPEND_SUMMARY_TO_DAT_FILE)
        {
            ofstream xpfile(SUMMARY_FILE_NAME.c_str(), ofstream::app);
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
        }

        delete model;
    }
    else
    {
        // LP relaxation infeasible
        delete lpr_model;
    }

    delete instance;
    return 0;
}
