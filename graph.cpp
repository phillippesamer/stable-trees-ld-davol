#include "graph.h"

Graph::Graph()
{
    using_matrix = false;
    using_lemon = false;
    num_vertices = 0;
    num_edges = 0;
}

Graph::Graph(long n, long m)
{
    using_matrix = false;
    using_lemon = false;
    num_vertices = n;
    num_edges = m;

    this->s.reserve(m);
    this->t.reserve(m);
    this->w.reserve(m);
}

Graph::~Graph()
{
    s.clear();
    t.clear();
    w.clear();

	if (using_matrix)
        free_index_matrix();

    if (using_lemon)
    {
        lemon_vertices.clear();
        lemon_edges.clear();
        delete lemon_weight;
        delete lemon_graph;
    }
}

void Graph::init_index_matrix()
{
	using_matrix = true;

    index_matrix = new long*[num_vertices];
    for (long i=0; i<num_vertices; ++i)
    {
        index_matrix[i] = new long[num_vertices];
        for (long j=0; j<num_vertices; ++j)
            index_matrix[i][j] = -1;
    }
}

void Graph::free_index_matrix()
{
    for (long i = 0; i<num_vertices; ++i)
        delete[] index_matrix[i];

    delete[] index_matrix;
}

void Graph::init_lemon()
{
    using_lemon = true;

    lemon_graph = new ListGraph();

    lemon_vertices.reserve(num_vertices);
    for (long i=0; i<num_vertices; ++i)
        lemon_vertices.push_back(lemon_graph->addNode());

    lemon_edges.reserve(num_edges);
    lemon_weight = new ListGraph::EdgeMap<long>(*lemon_graph);
}

void Graph::update_single_weight(long idx, long new_weight)
{
    // weight in the edge list
    this->w[idx] = new_weight;

    // weight in lemon's adjacency list
    ListGraph::Edge e = this->lemon_edges[idx];
    (*lemon_weight)[e] = new_weight;
}

void Graph::update_all_weights(vector<long> new_weights)
{
    // weight in the edge list
    this->w.clear();
    this->w = vector<long>(new_weights);

    // weight in lemon's adjacency list
    for (long i=0; i<num_edges; ++i)
    {
        ListGraph::Edge e = this->lemon_edges[i];
        (*lemon_weight)[e] = new_weights[i];
    }
}

bool Graph::mst()
{
    // 1. use lemon object to calculate mst via efficient implementation of
    // kruskal's algorithm

    // 2. store mst cost and solution in this object

    /*
    this->mst_cost = ...;

    this->mst_edges.clear();
    this->mst_edges = vector<long>();
    ...
    */
}

bool Graph::lemon_test_adj(ListGraph &g, ListGraph::Node &x, ListGraph::Node &y)
{
    // auxiliary function to test adjacency in the LEMON data structure
    for (ListGraph::IncEdgeIt e(g, x); e != INVALID; ++e)
    {
        //cout << "testing if " << g.id(y) << " is equal to " << g.id(g.v(e)) << " or " << g.id(g.u(e)) << endl;
        if ( g.id(g.v(e)) == g.id(y) || g.id(g.u(e)) == g.id(y))
            return true;
    }
    return false;
}
