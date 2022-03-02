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

// parameters controlled by the user
class UFL_parms {
public:
  std::string fdata; // file with the data
  std::string dualfile; // file with an initial dual solution
  std::string dual_savefile; // file to save final dual solution
  std::string int_savefile; // file to save primal integer solution
  int h_iter; // number of times that the primal heuristic will be run
  // after termination of the volume algorithm
   
  UFL_parms(const char* filename);
  ~UFL_parms() {}
};


/***
 * \file lddda_vol.h
 * 
 * Module especializing the volume algorithm implementation available in the
 * COIN-OR Vol project (see https://github.com/coin-or/Vol). Also inherits
 * from the LDDA class which approximates the lagrangean decomposition bound
 * with a dual ascent approach.
 * 
 * NB: Vol is licensed under the EPL license.
 * 
 * \author Phillippe Samer <phillippes@gmail.com>
 * \date 03.02.2022
 */
class LDDAVolume : public VOL_user_hooks, public LDDA
{
public:
    LDDAVolume(IO*, KStabModel*);
    LDDAVolume(IO*, KStabModel*, vector<long>);
    virtual ~LDDAVolume();
    LDDAVolume(const LDDAVolume &);

    bool run_volume();
    void UFL_read_data(const char*);

  // for all hooks: return value of -1 means that volume should quit

  // compute reduced costs
   int compute_rc(const VOL_dvector& u, VOL_dvector& rc);

  // solve relaxed problem
   int solve_subproblem(const VOL_dvector& u, const VOL_dvector& rc,
			double& lcost, VOL_dvector& x, VOL_dvector&v,
			double& pcost);

  // primal heuristic
   int heuristics(const VOL_problem& p, 
		  const VOL_dvector& x, double& heur_val);

  // original data for uncapacitated facility location
public: 
  VOL_dvector fcost; // cost for opening facilities
  VOL_dvector dist; // cost for connecting a customer to a facility
  VOL_dvector fix; // vector saying if some variables should be fixed
  // if fix=-1 nothing is fixed
  int ncust, nloc; // number of customers, number of locations
  VOL_ivector ix;   // best integer feasible solution so far
  double      icost;  // value of best integer feasible solution 
};

#endif
