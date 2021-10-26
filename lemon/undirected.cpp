#include <iostream>
#include <vector>
#include <lemon/list_graph.h>
#include <lemon/kruskal.h>
#include <lemon/time_measure.h>
#include <lemon/core.h>

using namespace lemon;
using namespace std;

bool lemon_test_adj(ListGraph &g, ListGraph::Node &x, ListGraph::Node &y)
{
    for (ListGraph::IncEdgeIt e(g, x); e != INVALID; ++e)
    {
        //cout << "testing if " << g.id(y) << " is equal to " << g.id(g.v(e)) << " or " << g.id(g.u(e)) << endl;
        if ( g.id(g.v(e)) == g.id(y) || g.id(g.u(e)) == g.id(y))
            return true;
    }
    return false;
}

int main()
{
    long num_vertices = 5;
    long big_M = 100000;

    ListGraph l_graph;
    vector<ListGraph::Node> l_vertices;
    vector<ListGraph::Edge> l_edges;
    ListGraph::EdgeMap<long> l_weight(l_graph);

    for (long i=0; i<num_vertices; ++i)
        l_vertices.push_back(l_graph.addNode());

    for (long i=0; i<num_vertices; ++i)
    {
        ListGraph::Edge e = l_graph.addEdge(l_vertices[i], l_vertices[i+1 % 9]);
        l_edges.push_back(e);
        l_weight[e] = num_vertices*i;
    }

    // adding all other edges with much higher weights
    for (ListGraph::NodeIt u(l_graph); u != INVALID; ++u)
    {
        for (ListGraph::NodeIt v(l_graph); v != INVALID; ++v)
        {
            if ( u!=v && !lemon_test_adj(l_graph, u, v) )
            {
                ListGraph::Edge e = l_graph.addEdge(u,v);
                l_edges.push_back(e);
                l_weight[e] = big_M;
            }
        }
    }

    cout << "n = " << countNodes(l_graph) << endl;
    cout << "m = " << countEdges(l_graph) << endl;

    // iterating over edges
    for (ListGraph::EdgeIt e_it(l_graph); e_it != INVALID; ++e_it)
    {
        cout << "edge " << l_graph.id(e_it) << " is {" << l_graph.id(l_graph.u(e_it)) << "," << l_graph.id(l_graph.v(e_it)) << "}, ";
        cout << "weight = " << l_weight[e_it] << endl;
    }

    // iterating over vertices and their neighbourhood
    for (ListGraph::NodeIt vertex(l_graph); vertex != INVALID; ++vertex)
    {
        int cnt = 0;
        for (ListGraph::IncEdgeIt e_it(l_graph,vertex); e_it != INVALID; ++e_it)
            cnt++;

        cout << "deg(" << l_graph.id(vertex) << ") = " << cnt << endl;
    }

    // mst algorithm
    Timer mst_timer;
    mst_timer.start();
    vector<ListGraph::Edge> tree_vector;
    long mst_weight = kruskal(l_graph, l_weight, back_inserter(tree_vector));
    mst_timer.halt();

    cout << "MST weight = " << mst_weight << endl;
    cout << "MST edges: " << endl;
    for (unsigned i = 0; i != tree_vector.size(); ++i)
        cout << "{" << l_graph.id(l_graph.u(tree_vector[i])) << ", " << l_graph.id(l_graph.v(tree_vector[i])) << "}" << endl;
    cout << "Elapsed time: " << mst_timer.realTime() << endl;

    /////////////////////////////////////////////////////////////////////////////////
    // copying a graph
    ListGraph l_graph_2;
    GraphCopy<ListGraph, ListGraph> cg(l_graph, l_graph_2);
    // save vertex correspondence: must be from old to new
    ListGraph::NodeMap<ListGraph::Node> nr(l_graph);
    cg.nodeRef(nr);
    // save inverted edge correspondences: must be from new to old
    ListGraph::EdgeMap<ListGraph::Edge> ecr(l_graph_2);
    cg.edgeCrossRef(ecr);
    // copy edge map
    ListGraph::EdgeMap<long> l_weight_2(l_graph_2);
    cg.edgeMap(l_weight, l_weight_2);
    // Execute copying
    cg.run();


    // updating weights
    /*
    for (ListGraph::EdgeIt e(l_graph_2); e != INVALID; ++e)
        if (l_weight_2[e] < big_M)
            l_weight_2[e] *= 10;
    */

    // removing expensive edges
    /*
    for (ListGraph::EdgeIt e(l_graph_2); e != INVALID; ++e)
        if (l_weight_2[e] >= big_M)
            l_graph_2.erase(e);
    */

    cout << "contracting edge {" << l_graph_2.id(nr[l_vertices[3]]) << "," << l_graph_2.id(nr[l_vertices[2]]) << "}" << endl;
    l_graph_2.contract(nr[l_vertices[3]], nr[l_vertices[2]], false);


    cout << "n_2 = " << countNodes(l_graph_2) << endl;
    cout << "m_2 = " << countEdges(l_graph_2) << endl;

    // iterating over edges
    for (ListGraph::EdgeIt e_it(l_graph_2); e_it != INVALID; ++e_it)
    {
        cout << "edge " << l_graph_2.id(e_it) << " is {" << l_graph_2.id(l_graph_2.u(e_it)) << "," << l_graph_2.id(l_graph_2.v(e_it)) << "}, ";
        cout << "weight = " << l_weight_2[e_it] << ", ";
        cout << "aka {" << l_graph.id(l_graph.u(ecr[e_it])) << ", " << l_graph.id(l_graph.v(ecr[e_it])) << "} in the original graph" << endl;
    }

    // iterating over vertices and their neighbourhood
    for (ListGraph::NodeIt vertex(l_graph_2); vertex != INVALID; ++vertex)
    {
        int cnt = 0;
        for (ListGraph::IncEdgeIt e_it(l_graph_2,vertex); e_it != INVALID; ++e_it)
            cnt++;

        cout << "deg(" << l_graph_2.id(vertex) << ") = " << cnt << endl;
    }

    // mst algorithm on graph 2
    mst_timer.reset();
    mst_timer.start();
    vector<ListGraph::Edge> tree_vector_2;
    long mst_weight_2 = kruskal(l_graph_2, l_weight_2, back_inserter(tree_vector_2));
    mst_timer.halt();

    cout << "MST weight = " << mst_weight_2 << endl;
    cout << "MST edges: " << endl;
    for (unsigned i = 0; i != tree_vector_2.size(); ++i)
    {
        cout << "{" << l_graph_2.id(l_graph_2.u(tree_vector_2[i])) << ", " << l_graph_2.id(l_graph_2.v(tree_vector_2[i])) << "} ";
        cout << "aka {" << l_graph.id(l_graph.u(ecr[tree_vector_2[i]])) << ", " << l_graph.id(l_graph.v(ecr[tree_vector_2[i]])) << "} in the original graph" << endl;

    }
    cout << "Elapsed time: " << mst_timer.realTime() << endl;


    // mst algorithm on graph 1
    cout << "n = " << countNodes(l_graph) << endl;
    cout << "m = " << countEdges(l_graph) << endl;

    mst_timer.reset();
    mst_timer.start();
    tree_vector.clear();
    mst_weight = kruskal(l_graph, l_weight, back_inserter(tree_vector));
    mst_timer.halt();

    cout << "MST weight = " << mst_weight << endl;
    cout << "MST edges: " << endl;
    for (unsigned i = 0; i != tree_vector.size(); ++i)
        cout << "{" << l_graph.id(l_graph.u(tree_vector[i])) << ", " << l_graph.id(l_graph.v(tree_vector[i])) << "}" << endl;
    cout << "Elapsed time: " << mst_timer.realTime() << endl;
}
