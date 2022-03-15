/***
 * \file main.cpp
 * 
 * Lagrangean Decomposition based dual ascent + volume algorithm to compute
 * lower bounds for minimum weight stable (or conflict-free) spanning trees.
 * 
 * We use the volume algorithm implementation available in the COIN-OR
 * Vol project (see https://github.com/coin-or/Vol).
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 02.03.2022
 */

#include "io.h"
#include "kstab_model.h"
#include "mstcc_model.h"
#include "ldda.h"
#include "ldda_vol.h"

#include <cstdlib>
#include <fstream>

using namespace std;

// execution switches
bool RUN_LDDA = true;
bool RUN_VOL = true;
bool INITIALIZE_VOL_FROM_LDDA = true;

double RUN_KSTAB_WITH_TIME_LIMIT = 3600;

bool RUN_STEEPEST_ASCENT_LDDA = false;
bool WRITE_LDDA_LOG_FILE = true;

bool APPEND_SUMMARY_TO_DAT_FILE = true;
string SUMMARY_FILE_NAME = string("xp4table.dat");

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cout << "usage: \t" << argv[0] << " [input file]" << endl << endl;
        return 0;
    }

    IO* instance = new IO();
    
    // parse given input file and look for errors in it
    if (instance->parse_input_file(string(argv[1])) == false)
    {
        cout << "unable to parse input file" << endl;
        delete instance;
        return 0;
    }

    // volume extension: primal bound from a max-weight mst
    instance->run_maxst();

    // 1. LP RELAXATION BOUND OF THE MSTCC NATURAL IP FORMULATION

    StableSpanningTreeModel *lpr_model = new StableSpanningTreeModel(instance);

    if ( !lpr_model->solve_lp_relax(false) )
    {
        // lp relaxation infeasible
        delete lpr_model;
    }
    else
    {
        // 2. TRIVIAL COMBINATORIAL BOUNDS FROM MIN WEIGHT SPANNING TREE AND KSTAB

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

        if (model->solution_status != AT_OPTIMUM)
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
        else
        {
            // lagrangean decomposition bound, approximated by dual ascent and/or volume
            LDDAVolume *lagrangean = new LDDAVolume(instance, model);

            if (RUN_LDDA)
            {
                // 3. LAGRANGEAN DECOMPOSITION BOUND: APPROXIMATION BY DUAL ASCENT

                bool ldda_complete = lagrangean->dual_ascent(RUN_STEEPEST_ASCENT_LDDA);

                cout << endl << "ldda bound: ";
                if (ldda_complete)
                    cout << lagrangean->bound_log.back();
                else if (!ldda_complete && lagrangean->problem_solved)   // infeasible
                    cout << " x ";
                else
                    cout << " - ";
                cout << " (runtime " << fixed << lagrangean->runtime << ")" << endl << endl;

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
                    else if (!ldda_complete && lagrangean->problem_solved)   // infeasible
                        table_row << setw(10) << " x ";
                    else
                        table_row << setw(10) << " - ";
                    table_row << setw(5) << "  &  ";
                    table_row << setw(10) << fixed << lagrangean->runtime;
                    table_row << setw(5) << endl;
                }

            } // run ldda condition

            if (RUN_VOL && !lagrangean->problem_solved)
            {
                // 4. LAGRANGEAN DECOMPOSITION BOUND: APPROXIMATION BY THE VOLUME ALGORITHM

                cout << "maxwst primal bound: " << instance->get_maxst_weight() << endl;

                if (RUN_LDDA && INITIALIZE_VOL_FROM_LDDA)
                    lagrangean->initialize_multipliers( lagrangean->multipliers_log.back() );

                bool volume_complete = lagrangean->run_volume();

                cout << endl << "volume bound: ";
                if (volume_complete)
                    cout << lagrangean->volume_bound;
                else if (!volume_complete && lagrangean->problem_solved)   // infeasible
                    cout << " x ";
                else
                    cout << " - ";
                cout << " (" << lagrangean->volume_iterations << " iterations, runtime " << fixed << lagrangean->volume_runtime << ")" << endl;

            } // run volume condition

            delete lagrangean;

        } // kstab optimal condition

        // finally: actually write summary output to file (facilitates experiments with many instances)
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

    } // lp relaxation condition

    delete instance;
    return 0;
}
