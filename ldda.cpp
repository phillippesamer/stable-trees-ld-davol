#include "ldda.h"

LDDA::LDDA(IO *instance, Model *model)
{
    this->instance = instance;
    this->model = model;

    // initialize multipliers at zero
    multipliers_log = vector< vector<double> >();
    multipliers_log.push_back( vector<double>(instance->graph->num_edges, 0.0) );
}

LDDA::LDDA(IO *instance, Model *model, vector<double> initial_multipliers)
{
    this->instance = instance;
    this->model = model;

    // use given values as the initial multipliers
    multipliers_log = vector< vector<double> >();
    multipliers_log.push_back(initial_multipliers);
}

LDDA::~LDDA()
{
    multipliers_log.clear();
}

double LDDA::edge_deletion_bound()
{
    // copy lemon graph

    // delete edge

    return 0;
}

double LDDA::edge_contraction_bound()
{
    // copy lemon graph

    // contract edge

    return 0;
}

double LDDA::vertex_deletion_bound()
{
    return 0;
}

double LDDA::vertex_fix_bound()
{
    return 0;
}

double LDDA::dual_ascent()
{
    // main loop

    instance->graph->mst();

    model->solve(true);

    /*
    long m = instance->graph->num_edges;

    print_edge_weights();

    for (long i=0; i<m; ++i)
        instance->graph->update_single_weight(i,10*i);
    print_edge_weights();

    instance->graph->update_all_weights( vector<long>(m, -1) );
    print_edge_weights();

    */

    return 0;
}


void LDDA::print_edge_weights()
{
    // print current weights (in MST subproblem)
    cout << endl << "edge list weights: ";
    for (long idx=0; idx < instance->graph->num_edges; ++idx)
        cout << this->instance->graph->w[idx] << " ";

    cout << endl << "lemonlist weights: ";
    for (long idx=0; idx < instance->graph->num_edges; ++idx)
        cout << (*instance->graph->lemon_weight)[ instance->graph->lemon_edges[idx] ] << " ";
    cout << endl << endl;
}
