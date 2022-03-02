#include <cstdio>
#include <set>
#include <list>
#include <vector>
#include <cmath>
#include <sys/times.h>

#include "ldda_vol.h"

//****** UFL_parms
// reading parameters specific to facility location
UFL_parms::UFL_parms(const char *filename) :
   fdata(""), 
   h_iter(10)
{
   char s[500];
   FILE * file = fopen(filename, "r");
   if (!file) { 
      printf("Failure to open ufl datafile: %s\n ", filename);
      abort();
   }
   
   while (fgets(s, 500, file)) {
      const int len = strlen(s) - 1;
      if (s[len] == '\n')
     s[len] = 0;
      std::string ss;
      ss = s;
      
      if (ss.find("fdata") == 0) {
     int j = ss.find("=");
     int j1 = ss.length() - j + 1;
     fdata = ss.substr(j+1, j1);
     
      } else if (ss.find("dualfile") == 0) {
     int j = ss.find("=");
     int j1 = ss.length() - j + 1;
     dualfile = ss.substr(j+1, j1);
     
      } else if (ss.find("dual_savefile") == 0) {
     int j = ss.find("=");
     int j1 = ss.length() - j + 1;
     dual_savefile = ss.substr(j+1, j1);

      } else if (ss.find("int_savefile") == 0) {
     int j = ss.find("=");
     int j1 = ss.length() - j + 1;
     int_savefile = ss.substr(j+1, j1);
     
      } else if (ss.find("h_iter") == 0) {
     int i = ss.find("=");  
     h_iter = atoi(&s[i+1]);
      } 
   }
   fclose(file);
}

LDDAVolume::LDDAVolume(IO *instance, KStabModel *model)
: LDDA(instance, model), icost(DBL_MAX)
{

}

LDDAVolume::LDDAVolume(IO *instance, KStabModel *model, vector<long> initial_mult)
: LDDA(instance, model, initial_mult), icost(DBL_MAX)
{

}

LDDAVolume::~LDDAVolume()
{

}


LDDAVolume::LDDAVolume(const LDDAVolume &source)
: LDDA(source.instance, source.model)
{
    // copy constructor to guarantee deep copy of data vectors

    ncust = source.ncust;
    nloc = source.nloc;
    icost = source.icost;

    fcost.allocate(nloc);
    fcost = source.fcost;

    dist.allocate(nloc*ncust);
    dist = source.dist;

    fix.allocate(nloc);
    fix = source.fix;

    ix.allocate(nloc + nloc*ncust);
    ix = source.ix;
}


bool LDDAVolume::run_volume()
{
    //////////////////////////////////////////////////////////////////////////
    // TO DO: reflect graph changes inside LDDA in IO, Graph, models, etc...?
    // TO DO: remove UFL specific code, implement new functions
    // TO DO: give maxweightST primal bound to volume
    // TO DO: set correct filename of the initial multipliers on .par file

    //////////////////////////////////////////////////////////////////////////

    // read in problem specific parameters and initialize data structures
    UFL_parms ufl_par("ldda_vol.par");   // reads parameter file
    this->icost = DBL_MAX;
    this->UFL_read_data(ufl_par.fdata.c_str());   // reads data file

   // initializes the VOL_problem object (also reads the parameter file)
   VOL_problem volp("ldda_vol.par");
   volp.psize = this->nloc + (this->nloc * this->ncust); // # primal vars
   volp.dsize = this->ncust;    // # dual vars

   bool ifdual = false;
   if (ufl_par.dualfile.length() > 0) {
     // read dual solution
      ifdual = true;
      VOL_dvector& dinit = volp.dsol;
      dinit.allocate(volp.dsize);
      // read from file
      FILE * file = fopen(ufl_par.dualfile.c_str(), "r");
      if (!file) {
        printf("Failure to open file: %s\n ", ufl_par.dualfile.c_str());
        abort();
      }
      const int dsize = volp.dsize;
      int idummy;
      for (int i = 0; i < dsize; ++i) {
        fscanf(file, "%i%lf", &idummy, &dinit[i]);
        cout << "read x[" << i << "] = " << dinit[i] << ", ";
      }
      cout << endl;
      fclose(file);
   }

   // start time measurement
   double t0;
   struct tms timearr; clock_t tres;
   tres = times(&timearr); 
   t0 = timearr.tms_utime; 

   // invoke volume algorithm
   if (volp.solve(*this, ifdual) < 0)
   {
        printf("solve failed...\n");
   }
   else
   {
      const int n = this->nloc;
      const int m = this->ncust;

      // recompute the violation of the fractional primal solution
      VOL_dvector v(volp.dsize);
      const VOL_dvector& psol = volp.psol;
      v = 1;
      int i,j,k=n;
      for (j = 0; j < n; ++j)
      {
        for (i = 0; i < m; ++i)
        {
            v[i] -= psol[k];
            ++k;
        }
      }

      double vc = 0.0;
      for (i = 0; i < m; ++i)
         vc += fabs(v[i]);
      vc /= m;
      printf(" Average violation of final solution: %f\n", vc);

      if (ufl_par.dual_savefile.length() > 0)
      {
        // save dual solution
        FILE* file = fopen(ufl_par.dual_savefile.c_str(), "w");
        const VOL_dvector& u = volp.dsol;
        int n = u.size();
        int i;
        for (i = 0; i < n; ++i)
        {
            fprintf(file, "%8i  %f\n", i+1, u[i]);
        }
        fclose(file);
      }

      // run a couple more heuristics
      /*
      double heur_val;
      for (i = 0; i < ufl_par.h_iter; ++i)
      {
        heur_val = DBL_MAX;
        this->heuristics(volp, psol, heur_val);
      }
      */

      // save integer solution
      if (ufl_par.int_savefile.length() > 0)
      {
        FILE* file = fopen(ufl_par.int_savefile.c_str(), "w");
        const VOL_ivector& x = this->ix;
        const int n = this->nloc;
        const int m = this->ncust;
        int i,j,k=n;
        fprintf(file, "Open locations\n");
        for (i = 0; i < n; ++i)
        {
            if ( x[i]==1 )
               fprintf(file, "%8i\n", i+1);
        }
        fprintf(file, "Assignment of customers\n");
        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < m; ++j)
            {
                if ( x[k]==1 ) 
                    fprintf(file, "customer %i  location %i\n", j+1, i+1);
                ++k;
           }
        }

        fclose(file);
      }

   }
   printf(" Best integer solution value: %f\n", this->icost);
   printf(" Lower bound: %f\n", volp.value);
   
   // end time measurement
   tres = times(&timearr);
   double t = (timearr.tms_utime-t0)/100.;
   printf(" Total Time: %f secs\n", t);

   return true;
}





//############################################################################

//
void LDDAVolume::UFL_read_data(const char* fname)
{

   FILE * file = fopen(fname, "r");
   if (!file) {
      printf("Failure to open ufl datafile: %s\n ", fname);
      abort();
   }


   VOL_dvector& fcost = this->fcost;
   VOL_dvector& dist = this->dist;

   int& nloc = this->nloc;
   int& ncust = this->ncust;
   int len;
#if 1
   char s[500];
   fgets(s, 500, file);
   len = strlen(s) - 1;
   if (s[len] == '\n')
      s[len] = 0;
   // read number of locations and number of customers
   sscanf(s,"%d%d",&nloc,&ncust);
 
   fcost.allocate(nloc);
   dist.allocate(nloc*ncust);
   double cost;
   int i,j,k;
   // read location costs
   for (i = 0; i < nloc; ++i) { 
     fgets(s, 500, file);
     len = strlen(s) - 1;
     if (s[len] == '\n')
	s[len] = 0;
     sscanf(s,"%lf",&cost);
     fcost[i]=cost;
   }
   dist=1.e7;
   while(fgets(s, 500, file)){
     len = strlen(s) - 1;
     if (s[len] == '\n')
	s[len] = 0;
     // read cost of serving a customer from a partucular location
     k=sscanf(s,"%d%d%lf",&i,&j,&cost);
     if(k!=3) break;
     if(i==-1)break;
     dist[(i-1)*ncust + j-1]=cost;
   }
#else
   fscanf(file, "%i%i", &ncust, &nloc);
   fcost.allocate(nloc);
   dist.allocate(nloc*ncust);
   int i,j;
   for ( j=0; j<ncust; ++j){
     for ( i=0; i<nloc; ++i){
       fscanf(file, "%f", &dist[i*ncust + j]);
     }
   }
   for ( i=0; i<nloc; ++i)
     fscanf(file, "%f", &fcost[i]);
#endif
   fclose(file);

   this->fix.allocate(nloc);
   this->fix = -1;
}

//############################################################################

//###### USER HOOKS
// compute reduced costs
int
LDDAVolume::compute_rc(const VOL_dvector& u, VOL_dvector& rc)
{
    /*** 
     * read current multipliers in vector u and determine the objective
     * coefficients of the next lagrangean subproblem (e.g. c_ij - u_j in UFL)
     * 
     * TO DO: set rc[i] = w_i - u[i] for i in 0...|E|-1
     * TO DO: set rc[i] =       u[i] for i in |E|...2*|E|-1
     */

   int i,j,k=0;
   for ( i=0; i < nloc; i++)
   {
     rc[i]=fcost[i];
     for (j = 0; j < ncust; ++j)
     {
       rc[nloc+k]= dist[k] - u[j];
       ++k;
     }
   }

   return 0;
}

// IN: dual vector u
//     reduced costs rc
// OUT: primal solution to the Lagrangian subproblem (x)
//      optimal value of Lagrangian subproblem (lcost)
//      v = difference between the rhs and lhs when substituting
//                  x into the relaxed constraints (v)
//      objective value of x substituted into the original problem (pcost)
// return value: -1 (volume should quit) 0 ow

int LDDAVolume::solve_subproblem( const VOL_dvector& u, 
                           const VOL_dvector& rc,
                           double& lcost,
                           VOL_dvector& x,
                           VOL_dvector& v,
                           double& pcost )
{
    /* TO DO: solve lagrangean subproblems with reduced costs rc
     * TO DO: save solution as dvector x
     * TO DO: determine lagrangean obj (lcost) and original obj (pcost)
     * TO DO: determine mismatch vector v =  "0 -(x-y)" with the subproblem solutions
     */

   int i,j;

   lcost = 0.0;
   for (i = 0; i < ncust; ++i) {
      lcost += u[i];
      v[i]=1;
   }

   VOL_ivector sol(nloc + nloc*ncust);

   // produce a primal solution of the relaxed problem
   const double * rdist = rc.v + nloc;
   double sum;
   int k=0, k1=0;
   double value=0.;
   int xi;
   for ( i=0; i < nloc; ++i ) {
     sum=0.;
     for ( j=0; j < ncust; ++j ) {
       if ( rdist[k]<0. ) sum+=rdist[k];
       ++k;
     }
     if (fix[i]==0) xi=0;
     else 
       if (fix[i]==1) xi=1;
       else 
	 if ( fcost[i]+sum >= 0. ) xi=0;
	 else xi=1;
     sol[i]=xi;
     value+=(fcost[i]+sum)*xi;
     for ( j=0; j < ncust; ++j ) {
       if ( rdist[k1] < 0. ) sol[nloc+k1]=xi;
       else sol[nloc+k1]=0;
       ++k1;
     }
   }

   lcost += value;

   pcost = 0.0;
   x = 0.0;
   for (i = 0; i < nloc; ++i) {
     pcost += fcost[i] * sol[i];
     x[i] = sol[i];
   }

   k = 0;
   for ( i=0; i < nloc; i++){
     for ( j=0; j < ncust; j++){
       x[nloc+k]=sol[nloc+k];
       pcost+= dist[k]*sol[nloc+k];
       v[j]-=sol[nloc+k];
       ++k;
     }
   }

   return 0;
}


// IN:  fractional primal solution (x),
//      best feasible soln value so far (icost)
// OUT: integral primal solution (this->ix) if better than best so far
//      and primal value (this->icost)
// returns -1 if Volume should stop, 0/1 if feasible solution wasn't/was
//         found.
//
int LDDAVolume::heuristics( const VOL_problem& p,
            		 const VOL_dvector& x,
                     double& new_icost)
{
    // TO DO: check if solution is integral, and then feasible
    // TO DO: (MAYBE) if solution is not integral, round some of the
    // mismatching vars and check feasibility

    return 0;
}
