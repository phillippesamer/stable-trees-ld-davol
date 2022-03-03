#ifndef __LDDAVOL_H__
#define __LDDAVOL_H__

#include "kstab_model.h"
#include "ldda.h"

#include <cstdio>
#include <cfloat>
#include <cstring>
#include <string>

// COIN-OR Vol header
#include "VolVolume.hpp"

/***
 * \file lddda_vol.h
 * 
 * Module especializing the volume algorithm implementation available in the
 * COIN-OR Vol project (see https://github.com/coin-or/Vol). Also inherits
 * from the LDDA class which approximates the lagrangean decomposition bound
 * with a dual ascent approach.
 * 
 * NB: COIN-OR Vol is licensed under the EPL license.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 03.02.2022
 */
class LDDAVolume : public VOL_user_hooks, public LDDA
{
public:
    LDDAVolume(IO*, KStabModel*);
    virtual ~LDDAVolume();

    bool read_volume_config_file(string);

    bool run_volume();
    double volume_runtime;
    long volume_iterations;
    vector<double> volume_bound_log;
    vector< vector<double> > volume_multipliers_log;  // lambda[iteration][var_index]

    // user hooks in COIN-OR Vol (NB! return -1 to quit volume)
    int compute_rc(const VOL_dvector&, VOL_dvector&);
    int solve_subproblem(const VOL_dvector&, const VOL_dvector&, double&, VOL_dvector&, VOL_dvector&, double&);
    int heuristics(const VOL_problem&, const VOL_dvector&, double&);

private:
    VOL_problem* volp;      // COIN-OR Vol problem description

    string config_file;    // config file
    bool initialized_mult; // true if the user defined an initial multipliers file

    string dualfile;       // initial multipliers
    string dual_savefile;  // file to save final multipliers
    long h_iter;            // heuristic runs after volume
};

#endif
