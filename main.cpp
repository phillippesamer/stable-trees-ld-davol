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
bool RUN_PREPROCESSING = true;
bool RUN_LDDA = true;
bool RUN_VOL = true;
bool INITIALIZE_VOL_FROM_LDDA = true;

double RUN_KSTAB_WITH_TIME_LIMIT = 1800;

bool RUN_STEEPEST_ASCENT_LDDA = false;
bool WRITE_LDDA_LOG_FILE = false;

bool WRITE_XP_FILE = true;
bool WRITE_XP_SECTION_SIMPLE_BOUNDS = true;
bool WRITE_XP_SECTION_LDDA_COLUMNS = true;
bool WRITE_XP_SECTION_VOL_COLUMNS = true;
bool WRITE_XP_VOL_TIME_INCLUDING_LDDA = true;
bool WRITE_XP_ROUNDING_BOUNDS_UP = true;
bool WRITE_XP_TAKING_MAX_OVER_KSTAB_MST_AND_LDDA_VOL = true;
string XP_FILE_NAME = string("xp12table.dat");

string trim_zeros(double, int);

int main(int argc, char **argv)
{
    if (argc < 2 || argc > 3)
    {
        cout << "usage: \t" << argv[0] << " [input file path] [optimum (optional)]" << endl << endl;
        cout << "the optimal value (when available) is only used to print how tight the Lagrangean bound is" << endl << endl;
        return 0;
    }

    IO* instance = new IO();
    bool instance_opt_given = (argc == 3);

    // 0. PARSE INPUT FILE AND RUN PREPROCESSING ALGORITHM
    if (instance->parse_input_file(string(argv[1]), RUN_PREPROCESSING) == false)
    {
        cout << "unable to parse input file" << endl;
        delete instance;
        return 0;
    }

    if (instance->problem_solved || instance->problem_infeasible)
    {
        if (instance->problem_solved)
        {
            instance->run_mst();
            cout << "mst weight: " << instance->get_mst_weight() + instance->objective_offset
                 << " (runtime " << fixed << instance->get_mst_runtime() << ")" << endl;
        }

        // output for xp file
        stringstream table_row;
        if (WRITE_XP_FILE)
        {
            table_row << left;
            table_row << setw(30) << instance->instance_id_trimmed;
            table_row << setw(5) << "  &  ";

            if (instance_opt_given)
                table_row << setw(10) << argv[2];
            else
                table_row << setw(10) << "";
            table_row << setw(6) << "  &&  ";

            if (WRITE_XP_SECTION_SIMPLE_BOUNDS)
            {
                table_row << setw(25) << "";
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << "";
                table_row << setw(6) << "  &&  ";
                if (instance->problem_solved)
                {
                    table_row << setw(10) << fixed << setprecision(0) << instance->get_mst_weight() + instance->objective_offset;
                    table_row << setw(5) << "  &  ";
                    table_row << setw(10) << fixed << setprecision(1) << instance->get_mst_runtime() + instance->preprocessing_runtime;
                }
                else
                {
                    table_row << setw(10) << "x";
                    table_row << setw(5) << "  &  ";
                    table_row << setw(10) << fixed << setprecision(1) << instance->preprocessing_runtime;
                }
                table_row << setw(6) << "  &&  ";
                table_row << setw(10) << "";
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << "";
                table_row << setw(6) << "  &&  ";
            }

            if (WRITE_XP_SECTION_LDDA_COLUMNS)
            {
                table_row << setw(10) << " - ";
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << " - ";
                table_row << setw(6) << "  &&  ";
            }

            if (WRITE_XP_SECTION_VOL_COLUMNS)
            {
                table_row << setw(10) << " - ";
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << " - ";
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << " - ";
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << "   \\\\  ";
            }

            table_row << setw(5) << "";

            ofstream xpfile(XP_FILE_NAME.c_str(), ofstream::app);
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

        delete instance;
        return 0;
    }

    // if running heuristics in volume: give primal bound from a max-weight mst
    // instance->run_maxst();

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

        double lp_bound_for_xp = WRITE_XP_ROUNDING_BOUNDS_UP ? 
            ceil(lpr_model->lp_bound + instance->objective_offset) :
            lpr_model->lp_bound + instance->objective_offset;
        double kstab_bound_for_xp = model->solution_weight + instance->objective_offset;
        double mst_bound_for_xp = instance->get_mst_weight() + instance->objective_offset;

        cout << "_____________________________________________________________________________" << endl << endl;

        cout << "kstab bound: " << model->solution_dualbound + instance->objective_offset
             << " (runtime " << fixed << model->solution_runtime << ")";
        if (model->solution_status != AT_OPTIMUM)
            cout << " *** NOT OPTIMAL ***";
        cout << endl;

        cout << "mst bound: " << instance->get_mst_weight() + instance->objective_offset
             << " (runtime " << fixed << instance->get_mst_runtime() << ")" << endl;

        cout << "lp bound: " << lpr_model->lp_bound + instance->objective_offset
             << " (" << lpr_model->lp_passes << " passes,"
             << " runtime " << fixed << lpr_model->lp_runtime << ")" << endl;
        cout << "_____________________________________________________________________________" << endl << endl;

        // summary output file (for experiments with many instances)
        stringstream table_row;
        if (WRITE_XP_FILE)
        {
            table_row << left;
            table_row << setw(30) << instance->instance_id_trimmed;
            table_row << setw(5) << "  &  ";

            if (instance_opt_given)
                table_row << setw(10) << argv[2];
            else
                table_row << setw(10) << "";
            table_row << setw(6) << "  &&  ";

            if (WRITE_XP_SECTION_SIMPLE_BOUNDS)
            {
                if (model->solution_status == AT_OPTIMUM)
                    table_row << setw(25) << model->solution_weight + instance->objective_offset;
                else if (model->solution_status == IS_INFEASIBLE)
                    table_row << setw(25) << "x";
                else
                {
                    stringstream tmp_str;
                    tmp_str << model->solution_weight + instance->objective_offset << " ?";
                    table_row << setw(25) << tmp_str.str();
                }
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << fixed << setprecision(1) << model->solution_runtime;
                table_row << setw(6) << "  &&  ";
                table_row << setw(10) << fixed << setprecision(0) << instance->get_mst_weight() + instance->objective_offset;
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << fixed << setprecision(1) << instance->get_mst_runtime();
                table_row << setw(6) << "  &&  ";
                table_row << setw(10) << trim_zeros(lp_bound_for_xp, 4);
                table_row << setw(5) << "  &  ";
                table_row << setw(10) << fixed << setprecision(1) << lpr_model->lp_runtime;
                table_row << setw(6) << "  &&  ";
            }
        }

        delete lpr_model;

        if (model->solution_status != AT_OPTIMUM)
        {
            /***
             * Either the kstab model could not be solved to optimality within
             * the time limit, or the preprocessing algorithm solved the problem
             * (proved infeasibility or removed all conflicts)
             */

            if (WRITE_XP_FILE)
            {
                if (WRITE_XP_SECTION_LDDA_COLUMNS)
                {
                    table_row << setw(10) << " - ";
                    table_row << setw(5) << "  &  ";
                    table_row << setw(10) << " - ";
                    table_row << setw(10) << " - ";
                    table_row << setw(5) << "  &&  ";
                }

                if (WRITE_XP_SECTION_VOL_COLUMNS)
                {
                    table_row << setw(10) << " - ";
                    table_row << setw(5) << "  &  ";
                    table_row << setw(10) << " - ";
                    table_row << setw(5) << "  &  ";
                    table_row << setw(10) << " - ";
                    table_row << setw(5) << "  &  ";
                    table_row << setw(10) << "   \\\\  ";
                }

                table_row << setw(5) << "";
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
                    cout << lagrangean->bound_log.back() + instance->objective_offset;
                else if (!ldda_complete && lagrangean->problem_solved)   // infeasible
                    cout << " x ";
                else
                    cout << " - ";
                cout << " (runtime " << fixed << lagrangean->runtime << ")" << endl << endl;

                double ldda_bound_for_xp = WRITE_XP_ROUNDING_BOUNDS_UP ?
                    ceil(lagrangean->bound_log.back() + instance->objective_offset) : 
                    lagrangean->bound_log.back() + instance->objective_offset;

                if (WRITE_XP_TAKING_MAX_OVER_KSTAB_MST_AND_LDDA_VOL)
                {
                    // printing the strongest between the kstab, mst, and ldda bounds
                    ldda_bound_for_xp = max(ldda_bound_for_xp, kstab_bound_for_xp);
                    ldda_bound_for_xp = max(ldda_bound_for_xp, mst_bound_for_xp);
                }

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

                if (WRITE_XP_FILE && WRITE_XP_SECTION_LDDA_COLUMNS)
                {
                    if (ldda_complete)
                        table_row << setw(10) << trim_zeros(ldda_bound_for_xp, 4);
                    else if (!ldda_complete && lagrangean->problem_solved)   // infeasible
                        table_row << setw(10) << " x ";
                    else
                        table_row << setw(10) << " - ";
                    table_row << setw(5) << "  &  ";
                    table_row << setw(10) << fixed << setprecision(1) << lagrangean->runtime;
                    table_row << setw(5) << "  &  ";
                    table_row << setw(10) << trim_zeros(100*(ldda_bound_for_xp - lp_bound_for_xp) / lp_bound_for_xp , 2);
                    table_row << setw(5) << "  &&  ";
                }

            } // run ldda condition

            if (RUN_VOL)
            {
                if (lagrangean->problem_solved)
                {
                    // no need to run the volume algorithm - dual ascent proved infeasibility
                    if (WRITE_XP_FILE && WRITE_XP_SECTION_VOL_COLUMNS)
                    {
                        table_row << setw(10) << " x ";
                        table_row << setw(5) << "  &  ";
                        table_row << setw(10) << fixed << setprecision(1) << lagrangean->runtime;
                        table_row << setw(5) << "  &  ";
                        table_row << setw(10) << " x ";
                        table_row << setw(5) << "  &  ";
                        table_row << setw(10) << "   \\\\  ";
                        table_row << setw(5) << "";
                    }
                }
                else
                {
                    // 4. LAGRANGEAN DECOMPOSITION BOUND: APPROXIMATION BY THE VOLUME ALGORITHM

                    // cout << "maxwst primal bound: " << instance->get_maxst_weight() << endl;

                    if (RUN_LDDA && INITIALIZE_VOL_FROM_LDDA)
                        lagrangean->initialize_multipliers( lagrangean->multipliers_log.back() );

                    bool volume_complete = lagrangean->run_volume();

                    cout << endl << "volume bound: ";
                    if (volume_complete)
                        cout << lagrangean->volume_bound + instance->objective_offset;
                    else if (!volume_complete && lagrangean->problem_solved)   // infeasible
                        cout << " x ";
                    else
                        cout << " - ";
                    cout << " (" << lagrangean->volume_iterations << " iterations, runtime " << fixed << lagrangean->volume_runtime << ")" << endl << endl;

                    double vol_bound_for_xp = WRITE_XP_ROUNDING_BOUNDS_UP ?
                        ceil(lagrangean->volume_bound + instance->objective_offset) : 
                        lagrangean->volume_bound + instance->objective_offset;

                    if (WRITE_XP_TAKING_MAX_OVER_KSTAB_MST_AND_LDDA_VOL)
                    {
                        // printing the strongest between the kstab, mst, and vol bounds
                        vol_bound_for_xp = max(vol_bound_for_xp, kstab_bound_for_xp);
                        vol_bound_for_xp = max(vol_bound_for_xp, mst_bound_for_xp);
                    }

                    if (WRITE_XP_FILE && WRITE_XP_SECTION_VOL_COLUMNS)
                    {
                        if (volume_complete)
                            table_row << setw(10) << trim_zeros(vol_bound_for_xp, 4);
                        else if (!volume_complete && lagrangean->problem_solved)   // infeasible
                            table_row << setw(10) << " x ";
                        else
                            table_row << setw(10) << " - ";
                        table_row << setw(5) << "  &  ";
                        double ldda_vol_runtime_on_xp = (WRITE_XP_VOL_TIME_INCLUDING_LDDA) ? (lagrangean->runtime + lagrangean->volume_runtime) : lagrangean->volume_runtime;
                        table_row << setw(10) << fixed << setprecision(1) << ldda_vol_runtime_on_xp;
                        table_row << setw(5) << "  &  ";
                        table_row << setw(10) << trim_zeros(100*(vol_bound_for_xp - lp_bound_for_xp) / lp_bound_for_xp , 2);
                        table_row << setw(5) << "  &  ";

                        if (instance_opt_given)
                        {
                            double opt = atof(argv[2]);
                            table_row << setw(10) << trim_zeros(100*(opt - vol_bound_for_xp) / opt , 1);
                        }

                        table_row << setw(10) << "   \\\\  ";
                        table_row << setw(5) << "";
                    }
                }

            } // run volume condition

            delete lagrangean;

        } // kstab optimal condition

        // finally: actually write summary output to file (facilitates experiments with many instances)
        if (WRITE_XP_FILE)
        {
            ofstream xpfile(XP_FILE_NAME.c_str(), ofstream::app);
            if (xpfile.is_open())
            {
                xpfile << table_row.str();
                xpfile << endl;
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

string trim_zeros(double value, int p)
{
    stringstream ss;
    ss << fixed << setprecision(p) << value;
    string str = ss.str();

    if(str.find('.') != string::npos)
    {
        // remove trailing zeros
        str = str.substr(0, str.find_last_not_of('0')+1);

        // remove the decimal point if it is now the last character
        if(str.find('.') == str.size()-1)
            str = str.substr(0, str.size()-1);
    }
    return str;
}
