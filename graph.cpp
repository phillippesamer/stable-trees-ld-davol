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
