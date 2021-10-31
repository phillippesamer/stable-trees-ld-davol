#ifndef _IO_H_
#define _IO_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <algorithm>

#include "graph.h"

using namespace std;

/// input and output functionality; includes a Graph object for main data structures
class IO
{
public:
    IO();
    ~IO();
    
    bool parse_gcclib(const char *);
    
private:
    friend class Model;

    // instance data
    long num_conflicts;
    string instance_id;

    Graph *graph;  // object with different graph representations

    vector< pair<long,long> > conflicts;  // conflicting edges (indexes only)
};

#endif
