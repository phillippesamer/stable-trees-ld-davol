#include <iostream>
#include <lemon/list_graph.h>

using namespace lemon;
using namespace std;

int main()
{
    ListDigraph g;
    ListDigraph::Node u = g.addNode();
    ListDigraph::Node v = g.addNode();
    ListDigraph::Arc  a = g.addArc(u, v);

    cout << "Hello World! This is LEMON library here." << endl;
    cout << "We have a directed graph with " << countNodes(g) << " nodes "
         << "and " << countArcs(g) << " arc." << endl;


    cout << "Arc " << g.id(a) << " is (" << g.id(g.source(a)) << "," << g.id(g.target(a)) << ")" << endl;

    for (ListDigraph::NodeIt n(g); n != INVALID; ++n)
        cout << g.id(n) << endl;
    
    for (ListDigraph::ArcIt a_it(g); a_it != INVALID; ++a_it)
            cout << "Arc " << g.id(a_it) << " is (" << g.id(g.source(a_it))
                 << "," << g.id(g.target(a_it)) << ")" << endl;

}
