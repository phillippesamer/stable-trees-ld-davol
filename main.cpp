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
#include "ldda.h"

#include <cstdlib>
#include <fstream>
#include <sys/time.h>       // for 'gettimeofday()'

using namespace std;

static struct timeval *clock_start;

void start_timer()
{
    // current clock time
    clock_start = (struct timeval *) malloc(sizeof(struct timeval) );
    gettimeofday(clock_start, 0);
}

void get_timer()
{
    struct timeval *clock_stop;
    clock_stop = (struct timeval *) malloc( sizeof(struct timeval) );
    gettimeofday(clock_stop, 0);

    unsigned long clock_time = 1.e6 * (clock_stop->tv_sec - clock_start->tv_sec) +
                                      (clock_stop->tv_usec - clock_start->tv_usec);

    printf( "main() says: runtime in seconds\n%.4f\n", ((double)clock_time / (double)1.e6) );

    free(clock_stop);
}

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

    KStabModel *model = new KStabModel(instance);
    LDDA *lagrangean = new LDDA(instance, model);

    /*
    if (model->solve_full_lp_relax(true))
    {
        cout << "lp_bound = " << model->full_lp_bound << "(runtime: " << model->full_lp_runtime << ")" << endl;
        // TODO: check acyclic solution
    }
    */

    //*
    start_timer();
    lagrangean->dual_ascent(false);
    get_timer();

    // write log to file
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
    //*/

    // clean up
    free(clock_start);
    delete lagrangean;
    delete instance;
    delete model;

    return 0;
}
