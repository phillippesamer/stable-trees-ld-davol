#include "io.h"

IO::IO()
{
    this->problem_solved = false;
    this->problem_infeasible = false;
    this->objective_offset = 0;
    this->num_conflicts = 0;
    this->preprocessing_runtime = 0;

    // only for consistency (if IO is used without calling parse_input_file)
    this->graph = new Graph();

    this->preprocessing_clock_start = (struct timeval *) malloc(sizeof(struct timeval));
    this->preprocessing_clock_stop = (struct timeval *) malloc(sizeof(struct timeval));
}

IO::~IO()
{
    delete graph;
    conflicts.clear();
    conflict_graph_adj_list.clear();

    free(preprocessing_clock_start);
    free(preprocessing_clock_stop);
}

// parse instance file
bool IO::parse_input_file(string filename, bool preprocessing)
{
    long num_vertices, num_edges;

    // check if this is an instance in the original benchmark (gcc extension)
    size_t dot_pos = filename.find_last_of(".");
    bool gcc_format = filename.find("gcc", dot_pos+1) != string::npos;

    ifstream input_fh(filename);
    
    if (input_fh.is_open())
    {
        // in gcc files, the first lines might have comments, followed by an instance id
        if (gcc_format)
        {
            string line;
            
            // skip comment lines
            do
            {
                getline(input_fh, line);
            }
            while(line.find("#") != string::npos);   // there is a trail

            instance_id.assign(line);
        }
        else
            instance_id.assign(filename);

        // trimmed instance id: contents after last slash and before the last dot
        size_t last_slash_pos = filename.find_last_of("/\\");
        instance_id_trimmed = filename.substr(last_slash_pos+1, dot_pos-1 - last_slash_pos);

        // 3 lines for cardinalities of vertex, edge and conflict sets
        input_fh >> num_vertices;
        input_fh >> num_edges;
        input_fh >> this->num_conflicts;

        // initialize graph (own structures and lemon object)
        delete graph;
        this->graph = new Graph(num_vertices,num_edges);
        this->graph->init_index_matrix();

        // conflict graph structures
        this->conflicts.reserve(this->num_conflicts);
        this->conflict_graph_adj_list.reserve(num_edges);
        this->conflict_graph_adj_list.insert(this->conflict_graph_adj_list.begin(), 
            num_edges, list<long>() );   // adds num_edges copies of an empty list
        
        // m lines for edges in the instance graph
        for (long line_idx=0; line_idx<num_edges; ++line_idx)
        {
            // read edge terminal vertices and weight: i j w  (incidentally, i<j)
            long i, j;
            double w;
            
            input_fh >> i;
            graph->s.push_back(i);

            input_fh >> j;
            graph->t.push_back(j);
            
            input_fh >> w;
            graph->w.push_back(w);

            graph->adj_list[i].push_back(j);
            graph->adj_list[j].push_back(i);
            
            // should never happen
            if (graph->index_matrix[i][j] >= 0 ||
                graph->index_matrix[j][i] >= 0 )
            {
                cerr << "ERROR: repeated edge in input file line " << line_idx << endl << endl;
                return false;
            }
            
            // store index of current edge
            graph->index_matrix[i][j] = line_idx;
            graph->index_matrix[j][i] = line_idx;
        }

        /* original instances (gcc format) include a single line per conflict,
         * but carrabs instances (cms format) include two lines for each (both
         * e1 e2 and e2 e1)
         */

        long conflict_lines_to_read = gcc_format ? num_conflicts : 2*num_conflicts;

        // p lines for conflicting edge pairs in the instance
        for (long line_idx=0; line_idx<conflict_lines_to_read; ++line_idx)
        {
            /* read terminal vertices from each edge in the pair:
             * a b c d (incidentally, a<b and c<d)
             */
            long a, b, c, d;
            input_fh >> a;
            input_fh >> b;
            input_fh >> c;
            input_fh >> d;
            
            long e1_index = graph->index_matrix[a][b];
            long e2_index = graph->index_matrix[c][d];

            // should never happen
            if (e1_index<0 || e2_index<0)
            {
                cerr << "ERROR: edge not found in conflict line " << line_idx << endl << endl;
                return false;
            }
            
            bool repeated = false;
            list<long>::iterator neighbour_it = this->conflict_graph_adj_list[e1_index].begin();
            while ( !repeated && neighbour_it != this->conflict_graph_adj_list[e1_index].end() )
            {
                if (*neighbour_it == e2_index)
                    repeated = true;

                ++neighbour_it;
            }

            // store conflict
            if (!repeated)
            {
                pair<long,long> edges = (a<c) ? make_pair(e1_index, e2_index) : make_pair(e2_index, e1_index);

                this->conflicts.push_back(edges);

                // update adjacency lists
                this->conflict_graph_adj_list[e1_index].push_back(e2_index);
                this->conflict_graph_adj_list[e2_index].push_back(e1_index);
            }
        }

        if ( (unsigned) num_conflicts != conflicts.size())
            cout << "num_conflicts = " << num_conflicts << ", conflict.size() = " << conflicts.size() << endl;
        
        input_fh.close();
    }
    else
    {
        cerr << "ERROR: could not open file (might not exist)." << endl;
        return false;
    }

    if (preprocessing)
    {
        start_timer();
        this->preprocess();
        this->preprocessing_runtime = total_time();
    }

    // lemon data structure initialization
    graph->init_lemon();

    return true;
}

bool IO::test_stability(vector<bool> &point)
{
    /// tests if an incidence vector satisfies all conflict constraints
    for ( vector< pair<long,long> >::iterator it = this->conflicts.begin();
          it != this->conflicts.end(); ++it )
    {
        long e1 = (*it).first;
        long e2 = (*it).second;
        if (point[e1] && point[e2])
            return false;
    }

    return true;
}

bool IO::test_acyclic(vector<bool> &point)
{
    /// tests if an incidence vector induces an acyclic subgraph

    bool acyclic = true;

    vector<bool> check = vector<bool>(graph->num_vertices, false);
    long checked_count = 0;

    // the search starts only once, if the subgraph is connected/acyclic
    long source = 0;
    while (source < graph->num_vertices && acyclic && checked_count < graph->num_vertices)
    {
        #ifdef DEBUG_DFS
            cout << "source = " << source << endl;
        #endif

        if (!check[source])
        {
            #ifdef DEBUG_DFS
                cout << "dfs starting!" << endl << endl;
            #endif
            dfs_checking_acyclic(source, -1, check, checked_count, point, acyclic);
        }

        ++source;
    }

    #ifdef DEBUG_DFS
        cout << "acyclic = " << acyclic << endl;
        cout << "check vector: { ";
        for (long i=0; i<graph->num_vertices; ++i) cout << check[i] << " ";
        cout << "}" << endl;
        cout << "checked_count = " << checked_count << endl;
    #endif

    return acyclic;
}

bool IO::test_acyclic_kstab(vector<bool> &point)
{
    /***
     * Tests if an incidence vector of num_vertices - 1 edges induces an
     * acyclic/connected subgraph.
     */

    bool acyclic = true;

    vector<bool> check = vector<bool>(graph->num_vertices, false);
    long checked_count = 0;

    // the search starts only once, if the subgraph is connected/acyclic
    dfs_checking_acyclic(0, -1, check, checked_count, point, acyclic);

    #ifdef DEBUG_DFS
        cout << "acyclic = " << (checked_count == graph->num_vertices) << endl;
        cout << "check vector: { ";
        for (long i=0; i<graph->num_vertices; ++i) cout << check[i] << " ";
        cout << "}" << endl;
        cout << "checked_count = " << checked_count << endl;
    #endif

    return (checked_count == graph->num_vertices);
}

void IO::dfs_checking_acyclic(long u,
                              long parent,
                              vector<bool> &check,
                              long &checked_count,
                              vector<bool> &active_edges,
                              bool &acyclic)
{
    /// dfs to check if the subgraph induced by active edges is acyclic

    check[u] = true;
    ++checked_count;

    list<long>::iterator it = graph->adj_list[u].begin();
    while (it != graph->adj_list[u].end() && acyclic)
    {
        long v = (*it);
        
        // considering only active edges
        long idx = graph->index_matrix[u][v];

        if (active_edges[idx])
        {
            if (!check[v])
                dfs_checking_acyclic(v, u, check, checked_count, active_edges, acyclic);
            else if (v != parent)
            {
                // vertex (not the antecessor) already visited: back edge!
                acyclic = false;
            }
        }

        ++it;
    }
}

void IO::start_timer()
{
    gettimeofday(preprocessing_clock_start, 0);
}

double IO::total_time()
{
    gettimeofday(preprocessing_clock_stop, 0);

    unsigned long clock_time = 1.e6 * (preprocessing_clock_stop->tv_sec - preprocessing_clock_start->tv_sec) +
                                      (preprocessing_clock_stop->tv_usec - preprocessing_clock_start->tv_usec);

    return (double) clock_time / (double)1.e6 ;
}


///////////////////////////////////////////////////////////////////////////////
// NOTE TO SELF: I implemented the code below in 2013, for the work presented
// in the 2015 paper in Optimization Letters with Sebastian. I see that I still
// used break and continue statements at the time...

// auxiliary procedure implementing depth-first search to store bridges
void IO::dfs_for_bridges(long u,
                         long &counter,
                         long parent,
                         vector<long> &bridges_idx,
                         bool *bridges_map,
                         long *chk,
                         long *lowest,
                         vector< vector<long> > &graph)
{
    chk[u] = lowest[u] = ++counter;
    
    for (unsigned i=0; i<graph[u].size(); ++i)
    {
        long v = graph[u].at(i);
        if (chk[v] < 0)
        {
            dfs_for_bridges(v, counter, u, bridges_idx, bridges_map, chk, lowest, graph);
            
            // v may lead to ancestor vertex
            lowest[u]  = min(lowest[u], lowest[v]);
            
            // closed vertex v: check if (u,v) is cut edge
            if (lowest[v] > chk[u])
            {
                long edge = this->graph->index_matrix[u][v];
                bridges_idx.push_back(edge);
                bridges_map[edge] = true;
            }
        }
        else
        {
            // v already visited, but check lowest
            if (v != parent)
                lowest[u] = min(lowest[u], chk[v]);
        }
    }
}

// analogous to the above dfs, but observe set and active edge lists during traversal
void IO::dfs_for_chain(long u,
                       long &counter,
                       long parent,
                       vector<long> &bridges_idx,
                       long *chk,
                       long *lowest,
                       vector< vector<long> > &graph,
                       bool *set_edges,
                       bool *inactive)
{
    chk[u] = lowest[u] = ++counter;
    
    for (unsigned i=0; i<graph[u].size(); ++i)
    {
        long v = graph[u].at(i);
        
        // skip edge (u,v) if inactive (temporarily removed in reduction chain)
        long idx = this->graph->index_matrix[u][v];
        if (inactive[idx] == true)
            continue;
        
        if (chk[v] < 0)
        {
            dfs_for_chain(v, counter, u, bridges_idx, chk, lowest, graph, set_edges, inactive);
            
            // v may lead to ancestor vertex
            lowest[u]  = min(lowest[u], lowest[v]);
            
            // closed vertex v: check if (u,v) is a cut edge, and store it if
            // not previously set in reduction chain
            if (lowest[v] > chk[u] && !set_edges[idx])
                bridges_idx.push_back(idx);
        }
        else
        {
            // v already visited, but check lowest
            if (v != parent)
                lowest[u] = min(lowest[u], chk[v]);
        }
    }
}


// fix bridges in original graph and remove corresponding conflicting edges
long IO::preprocess_bridges()
{
    bool updated = false;
    long fixed = 0;
    
    // caller procedure checks feasibility
    if (graph->num_edges <= graph->num_vertices-1)
        return 0;
    
    // auxiliary structures: will be overwritten several times as graph changes
    bool *bridges_map = new bool[graph->num_edges];
    long bridges_map_size = graph->num_edges;

    // checks visited vertices (stores discovery time) in dfs
    long *chk = new long[graph->num_vertices];
    long chk_size = graph->num_vertices;

    // earliest vertex seen in dfs
    long *lowest = new long[graph->num_vertices];
    long lowest_size = graph->num_vertices;
    
    // repeat while any edge is fixed or removed
    do
    {
        updated = false;
        
        // 1. additional data structures used in preprocessing bridges
        
        // original graph adjacency list
        vector< vector<long> > original_list;
        for (long i=0; i<graph->num_vertices; ++i)
            original_list.push_back(vector<long>());
    
        for (long e=0; e<graph->num_edges; ++e)
        {
            long i = graph->s[e];
            long j = graph->t[e];
        
            original_list[i].push_back(j);
            original_list[j].push_back(i);
        }
    
        // conflict graph adjacency list
        vector< vector<long> > conflict_list;
    
        for (long i=0; i<graph->num_edges; ++i)
            conflict_list.push_back(vector<long>());

        for (long p=0; p<num_conflicts; ++p)
        {
            pair<long,long> edges = conflicts.at(p);
            long e1 = edges.first;
            long e2 = edges.second;
            
            conflict_list[e1].push_back(e2);
            conflict_list[e2].push_back(e1);
        }
    
        // 2. find bridges (if any) in G using depth-first search
    
        vector<long> bridges_idx;
        memset(bridges_map, false, sizeof(bool)*bridges_map_size);
        memset(chk, -1, sizeof(long)*chk_size);
        memset(lowest, -1, sizeof(long)*lowest_size);
    
        long components = 0;
        long counter = -1;
        
        // shall execute dfs only once if G is connected
        for (long u=0; u<graph->num_vertices; ++u)
        {
            if (chk[u] < 0)
            {
                ++components;
                dfs_for_bridges(u, counter, -1, bridges_idx, bridges_map, chk, lowest, original_list);
            }
        }
    
        if (components > 1)
        {
            // infeasible instance (graph not connected) - clean up and return -1
            fixed = -1;
            updated = false;
        }
        else while (bridges_idx.size() > 0)
        {
            // 3. contract each bridge and remove conflicting edges
            updated = true;
            
            // NOTE: REMOVE ME LATER
            //cout << instance_id << " \t";
            //cout << graph->num_edges << " edges, " << bridges_idx.size() << " bridges in original graph" << endl;
            //for (unsigned b=0; b<bridges_idx.size(); ++b)
            //    cout << "(" << graph->s[bridges_idx[b]] << ", " << graph->t[bridges_idx[b]] << ")  ";
            //cout << endl;
            
            long bridge = bridges_idx.back();
            long v1 = graph->s[bridge];
            long v2 = graph->t[bridge];
            
            if (v1>v2)
            {
                long aux = v1;
                v1 = v2;
                v2 = aux;
            }
            
            // 4. cost of bridge shall be added to objective function
            objective_offset += graph->w[bridge];
            
            // 5. instance new data
            --(graph->num_vertices);               // two cut vertices are merged
            
            vector<long> new_bridges_idx;  // new index of bridges found
            vector<long> new_s;            // terminal vertex 1
            vector<long> new_t;            // terminal vertex 2
            vector<long> new_w;            // edge weight

            long **new_index_matrix = new long*[graph->num_vertices];  // adjacency matrix
            for (long i=0; i<graph->num_vertices; ++i)
            {
                new_index_matrix[i] = new long[graph->num_vertices];
                for (long j=0; j<graph->num_vertices; ++j)
                    new_index_matrix[i][j] = -1;
            }
            
            vector< pair<long,long> > new_conflicts;    // conflicting edges indexes
            vector< vector<long> > new_conflict_list;  // conflict graph adj list
            
            // 5.1 process each edge in original graph
            //cout << "graph had " << graph->num_edges << " edges" << endl;
            
            long new_num_edges = 0;
            for (long e=0; e<graph->num_edges; ++e)
            {
                // skip contracted edge
                if (e == bridge)
                {
                    //cout << "skipping contracted edge: #" << bridge << "("<< v1 << "," << v2 << ")" << endl;
                    ++fixed;
                    continue;
                }
                
                // skip conflicting edges
                bool skip_flag = false;
                for (vector<long>::iterator it = conflict_list[bridge].begin();
                     it != conflict_list[bridge].end(); ++it)
                {
                    if (e == *it)
                    {
                        //cout << "skipping conflicted edge: #" << e << "(" << graph->s[e] << "," << graph->t[e] << ")" << endl;
                        ++fixed;
                        
                        skip_flag = true;
                        break;   // impossible to use 'continue' here due to inner 'for'
                    }
                }
                if (skip_flag == true)
                    continue;
                
                // include edge (possibly updating vertex indexes): as
                // (v1,v2) is contracted, we keep the vertex with
                // smaller index, and refer to other k > v1 as k-1
                long i = graph->s[e];
                long j = graph->t[e];
                
                long new_i = (i <  v2) ? i  :
                            (i == v2) ? v1 : i-1;
                long new_j = (j <  v2) ? j  :
                            (j == v2) ? v1 : j-1;
                
                new_s.push_back(new_i);
                new_t.push_back(new_j);
                new_w.push_back(graph->w[e]);
                
                new_index_matrix[new_i][new_j] = new_num_edges;
                new_index_matrix[new_j][new_i] = new_num_edges;
                
                // keep track of bridges indexes as well 
                if (bridges_map[e])
                    new_bridges_idx.push_back(new_num_edges);
                
                ++new_num_edges;
            }
            
            //cout << "new graph has " << new_num_edges << " edges. old bridges are now: ";
            //for (unsigned b=0; b<new_bridges_idx.size(); ++b)
            //    cout << "(" << new_s[new_bridges_idx[b]] << ", " << new_t[new_bridges_idx[b]] << ")  ";
            //cout << endl;
            
            // 5.2 process each conflict (update conflict graph)
            
            for (long i=0; i<new_num_edges; ++i)
                new_conflict_list.push_back(vector<long>());
            
            long new_num_conflicts = 0;
            for (long p=0; p<num_conflicts; ++p)
            {
                pair<long,long> edges = conflicts.at(p);
                long e1 = edges.first;
                long e2 = edges.second;
                
                // remove conflict if any of the edges was fixed
                if (e1 == bridge || e2 == bridge)
                    continue;
                
                bool skip_flag = false;
                for (vector<long>::iterator it = conflict_list[bridge].begin();
                     it != conflict_list[bridge].end(); ++it)
                {
                    if (e1 == *it || e2 == *it)
                    {
                        skip_flag = true;
                        break;   // impossible to use 'continue' here due to inner 'for'
                    }
                }
                if (skip_flag)
                    continue;
                
                // find new indexes of conflicting edges and include conflict
                long new_e1_i = (graph->s[e1] <  v2) ? graph->s[e1]  :
                               (graph->s[e1] == v2) ? v1     : graph->s[e1]-1;
                long new_e1_j = (graph->t[e1] <  v2) ? graph->t[e1]  :
                               (graph->t[e1] == v2) ? v1     : graph->t[e1]-1;
                long new_e1 = new_index_matrix[new_e1_i][new_e1_j];
                
                long new_e2_i = (graph->s[e2] <  v2) ? graph->s[e2]  :
                               (graph->s[e2] == v2) ? v1     : graph->s[e2]-1;
                long new_e2_j = (graph->t[e2] <  v2) ? graph->t[e2]  :
                               (graph->t[e2] == v2) ? v1     : graph->t[e2]-1;
                long new_e2 = new_index_matrix[new_e2_i][new_e2_j];
                
                // TODO: REMOVE ME!
                if (new_e2 < 0 || new_e2 < 0)
                    cout << endl << endl << "UNEXPECTED: got -1 while trying to find new indexes of conflicting edges" << endl << endl;
                
                // include conflict in pairs list and adjacency list
                new_conflicts.push_back(make_pair(new_e1, new_e2));
                
                new_conflict_list[new_e1].push_back(new_e2);
                new_conflict_list[new_e2].push_back(new_e1);
                
                ++new_num_conflicts;
            }
            
            // TODO: REMOVE ME
            //long conflicting_count = conflict_list[bridge].size();
            //cout << "new_m = " << new_num_edges << " == " << graph->num_edges - (1+conflicting_count) << endl;
            for (long i=0; i<new_num_edges; ++i)
            {
                if (new_s[i] < 0 || new_s[i] >= graph->num_vertices ||
                    new_t[i] < 0 || new_t[i] >= graph->num_vertices)
                    cout << "edge joining vertex out of range [0, " << graph->num_vertices << "]" << endl;
            }
            
            // update this instance data structures
            graph->num_edges = new_num_edges;
            num_conflicts = new_num_conflicts;
            
            graph->s.clear();
            graph->t.clear();
            graph->w.clear();
            graph->s.assign(new_s.begin(), new_s.end());
            graph->t.assign(new_t.begin(), new_t.end());
            graph->w.assign(new_w.begin(), new_w.end());
            
            graph->free_index_matrix();
            graph->index_matrix = new_index_matrix;
            
            conflicts.clear();
            conflicts.assign(new_conflicts.begin(), new_conflicts.end());
            
            conflict_list.clear();
            conflict_list.assign(new_conflict_list.begin(), new_conflict_list.end());
            
            // adjust indexes of bridges for updated graph
            memset(bridges_map, false, sizeof(bool)*bridges_map_size);
            for (unsigned i=0; i<new_bridges_idx.size(); ++i)
                bridges_map[new_bridges_idx.at(i)] = true;
            bridges_idx.clear();
            bridges_idx.assign(new_bridges_idx.begin(), new_bridges_idx.end());
        }
        
    }
    while (updated && graph->num_edges>graph->num_vertices-1);
    
    delete[] chk;
    delete[] lowest;
    delete[] bridges_map;
    
    return fixed;
}

// auxiliary procedure, removing edge and corresponding conflicts
void IO::remove_edge(long idx)
{
    // 1. remove from adjacency matrix and fix every next index
    long u = graph->s[idx];
    long v = graph->t[idx];
    graph->index_matrix[u][v] = graph->index_matrix[v][u] = -1;
    
    for (long next=idx+1; next<graph->num_edges; ++next)
    {
        u = graph->s[next];
        v = graph->t[next];
        graph->index_matrix[u][v] = graph->index_matrix[v][u] = next-1;
    }
    
    // 2. remove from edge list
    graph->s.erase(graph->s.begin()+idx);
    graph->t.erase(graph->t.begin()+idx);
    graph->w.erase(graph->w.begin()+idx);
    
    --(graph->num_edges);
    
    // 3. remove conflicts (new structure overwrites that in instance)
    vector< pair<long,long> > new_conflicts;
    
    long new_num_conflicts = 0;
    for (long p=0; p<num_conflicts; ++p)
    {
        pair<long,long> edges = conflicts.at(p);
        long e1 = edges.first;
        long e2 = edges.second;
        
        // remove conflict if any of the edges was fixed
        if (e1 == idx || e2 == idx)
            continue;
        else
        {
            // new indexes of conflicting edges
            if (e1 > idx)
                --e1;
            
            if (e2 > idx)
                --e2;
            
            // include conflict in pairs list
            new_conflicts.push_back(make_pair(e1, e2));
            ++new_num_conflicts;
        }
    }
    
    conflicts.clear();
    conflicts.assign(new_conflicts.begin(), new_conflicts.end());
    
    num_conflicts = new_num_conflicts;
}

// auxiliary procedure, including conflict between edges
void IO::add_conflict(long idx1, long idx2)
{
    ++num_conflicts;
    conflicts.push_back(make_pair(idx1, idx2));
}

// finds if "conflict closure" of fixing given edge would disconnect graph
long IO::preprocess_chain()
{
    // 1. additional data structures for for bridge finding dfs
    
    // original graph adjacency list
    vector< vector<long> > original_list;
    for (long i=0; i<graph->num_vertices; ++i)
        original_list.push_back(vector<long>());

    for (long e=0; e<graph->num_edges; ++e)
    {
        long i = graph->s[e];
        long j = graph->t[e];
    
        original_list[i].push_back(j);
        original_list[j].push_back(i);
    }
    
    // conflict graph adjacency list
    vector< vector<long> > conflict_list;

    for (long i=0; i<graph->num_edges; ++i)
        conflict_list.push_back(vector<long>());

    for (long p=0; p<num_conflicts; ++p)
    {
        pair<long,long> edges = conflicts.at(p);
        long e1 = edges.first;
        long e2 = edges.second;
        
        conflict_list[e1].push_back(e2);
        conflict_list[e2].push_back(e1);
    }
    
    // next are some auxiliary structures, which are overwritten several times

    // checks visited vertices (stores discovery time) in dfs
    long *chk = new long[graph->num_vertices];

    // earliest vertex seen in dfs
    long *lowest = new long[graph->num_vertices];
    
    // used for temporary selection or removal of edges
    bool *inactive = new bool[graph->num_edges];
    bool *set_edges = new bool[graph->num_edges];
    
    // TODO: remove me
    //cout << "\t\t original graph has " << graph->num_edges << " edges:  ";
    //for (long e=0; e<graph->num_edges; ++e)
    //    cout << "#" << e << "(" << graph->s[e] << "," << graph->t[e] << ")  ";
    //cout << endl;
    //cout << "\t\t conflict graph has " << num_conflicts << " conflicts:   ";
    //for (long p=0; p<num_conflicts; ++p)
    //{
    //    pair<long,long> edges = conflicts.at(p);
    //    cout << "#" << edges.first << " vs #" << edges.second << "   ";
    //}
    //cout << endl;
    
    
    // 2. check connectivity fixing edges in chain, starting at first_e
    for (long first_e=0; first_e<graph->num_edges; ++first_e)
    {
        // TODO: REMOVE ME
        //if (conflict_list[first_e].empty())
        //    cout <<  "\t\t skipping conflict-free vertex #" << first_e << endl;
        
        // isolated vertices in the conflict graph cannot yield any reduction
        if (conflict_list[first_e].empty())
            continue;
        
        memset(set_edges, false, sizeof(bool)*graph->num_edges);
        memset(inactive, false, sizeof(bool)*graph->num_edges);
        
        // TODO: remove me
        //cout << "\t\t try fixing edge #" << first_e << " (" << graph->s[first_e] << "," << graph->t[first_e] << ")" << endl;
        //for (unsigned adj=0; adj<conflict_list[first_e].size(); ++adj)
        //    cout << "\t\t desactivating edge #" << conflict_list[first_e].at(adj) << " (" << graph->s[conflict_list[first_e].at(adj)] << "," << graph->t[conflict_list[first_e].at(adj)] << ")" << endl;
        
        // 3. temporarily select first_e, and reject conflicting ones
        long count_set_edges = 1;
        set_edges[first_e] = true;
        for (unsigned adj=0; adj<conflict_list[first_e].size(); ++adj)
            inactive[conflict_list[first_e].at(adj)] = true;
        
        // 4. perform bridge finding dfs, disregarding inactive edges
        vector<long> bridges_idx;
        memset(chk, -1, sizeof(long)*graph->num_vertices);
        memset(lowest, -1, sizeof(long)*graph->num_vertices);
    
        long components = 0;
        long counter = -1;
        
        // shall execute dfs only once if G is connected
        for (long u=0; u<graph->num_vertices; ++u)
        {
            if (chk[u] < 0)
            {
                ++components;
                dfs_for_chain(u, counter, -1, bridges_idx, chk, lowest, original_list, set_edges, inactive);
            }
        }
    
        // 5. selection leads to a disconnected graph: just remove first_e, done!
        if (components > 1)
        {
            // TODO: remove me
            //cout << "\t\t graph would be disconnected! removing edge #" << first_e << endl;
            
            remove_edge(first_e);
            
            // TODO: remove me
            //cout << "\t\t now the graph has " << graph->num_edges << " edges:  ";
            //for (long e=0; e<graph->num_edges; ++e)
            //    cout << "#" << e << "(" << graph->s[e] << "," << graph->t[e] << ")  ";
            //cout << endl;
            //cout << "\t\t now the conflict graph has " << num_conflicts << " conflicts:   ";
            //for (long p=0; p<num_conflicts; ++p)
            //{
            //    pair<long,long> edges = conflicts.at(p);
            //    cout << "#" << edges.first << " vs #" << edges.second << "   ";
            //}
            //cout << endl;
            
            delete[] chk;
            delete[] lowest;
            delete[] inactive;
            delete[] set_edges;
            
            return 1;
        }
        
        // 6. selection imply bridges in the original graph, which shall be fixed
        while (bridges_idx.size() > 0)
        {
            // NOTE: would there be a preferable choice among bridges?
            long bridge = bridges_idx.back();
            
            //cout << "\t\t\t\t implied fixing edge #" << bridge << " (" << graph->s[bridge] << "," << graph->t[bridge] << ")" << endl;
            //for (unsigned adj=0; adj<conflict_list[bridge].size(); ++adj)
            //    cout << "\t\t\t\t desactivating edge #" << conflict_list[bridge].at(adj) << " (" << graph->s[conflict_list[bridge].at(adj)] << "," << graph->t[conflict_list[bridge].at(adj)] << ")" << endl;
            
            // TODO: remove me
            if (inactive[bridge])
                cout << endl << endl << "UNEXPECTED: selected bridge was inactive!" << endl << endl;
            if (set_edges[bridge])
                cout << endl << endl << "UNEXPECTED: selected bridge was already set!" << endl << endl;
            
            // 7. temporary selection (conflicting edges could be inactive - ok!)
            ++count_set_edges;
            set_edges[bridge] = true;
            for (unsigned adj=0; adj<conflict_list[bridge].size(); ++adj)
                inactive[conflict_list[bridge].at(adj)] = true;
            
            // 8. perform a new dfs, finding bridges and connectivity status
            bridges_idx.clear();
            memset(chk, -1, sizeof(long)*graph->num_vertices);
            memset(lowest, -1, sizeof(long)*graph->num_vertices);
            
            components = 0;
            counter = -1;
            
            for (long u=0; u<graph->num_vertices; ++u)
            {
                if (chk[u] < 0)
                {
                    ++components;
                    dfs_for_chain(u, counter, -1, bridges_idx, chk, lowest, original_list, set_edges, inactive);
                }
            }
            
            // 9. selection leads to a disconnected graph: just remove first_e, done!
            if (components > 1)
            {
                // TODO: remove me
                //cout << "\t\t\t\t graph would be disconnected! removing edge #" << first_e << endl;
                
                remove_edge(first_e);
                
                // TODO: remove me
                //cout << "\t\t\t\t now the graph has " << graph->num_edges << " edges:  ";
                //for (long e=0; e<graph->num_edges; ++e)
                //    cout << "#" << e << "(" << graph->s[e] << "," << graph->t[e] << ")  ";
                //cout << endl;
                //cout << "\t\t\t\t now the conflict graph has " << num_conflicts << " conflicts:   ";
                //for (long p=0; p<num_conflicts; ++p)
                //{
                //    pair<long,long> edges = conflicts.at(p);
                //    cout << "#" << edges.first << " vs #" << edges.second << "   ";
                //}
                //cout << endl;
                
                delete[] chk;
                delete[] lowest;
                delete[] inactive;
                delete[] set_edges;

                return 1;
            }
            
            // 10. is current graph a (spanning, as it is connected) tree? stop chain
            if (count_set_edges == graph->num_vertices-1)
            {
                // stop chain as no more bridges will be checked
                bridges_idx.clear();
                
                // check if conflict constraints are satisfied
                bool feasible = true;
                for (long p=0; p<num_conflicts; ++p)
                {
                    pair<long,long> edges = conflicts.at(p);
                    long e1 = edges.first;
                    long e2 = edges.second;

                    // remove conflict if any of the edges was fixed
                    if (set_edges[e1] && set_edges[e2])
                    {
                        feasible = false;
                        break;
                    }
                }
                
                // store if integer feasible solution to use as a primal bound
                if (feasible)
                {
                    // TO DO: determine objective value and save (if its the new incumbent)
                    //cout << "[chain_reduction] found (primal) integer feasible solution " << endl;
                }
                else // no feasible spanning tree including this edge .: remove it
                {
                    // TO DO: remove me
                    //cout << "\t\t\t\t unique spanning tree is not feasible! removing edge #" << first_e << endl;
                    //cout << "[chain_reduction] removing edge because unique spanning tree is not feasible!" << endl;

                    remove_edge(first_e);

                    // TO DO: remove me
                    //cout << "\t\t\t\t now the graph has " << graph->num_edges << " edges:  ";
                    //for (long e=0; e<graph->num_edges; ++e)
                    //    cout << "#" << e << "(" << graph->s[e] << "," << graph->t[e] << ")  ";
                    //cout << endl;
                    //cout << "\t\t\t\t now the conflict graph has " << num_conflicts << " conflicts:   ";
                    //for (long p=0; p<num_conflicts; ++p)
                    //{
                    //    pair<long,long> edges = conflicts.at(p);
                    //    cout << "#" << edges.first << " vs #" << edges.second << "   ";
                    //}
                    //cout << endl;

                    delete[] chk;
                    delete[] lowest;
                    delete[] inactive;
                    delete[] set_edges;

                    return 1;
                }
            }
            
        } // fix another bridge (bridge_idx was updated properly)
        
        // if a cycle is implied, and first_e makes the problem infeasible .: just remove it
        /*
        if (count_set_edges >= 3)
        {
            // 11. perform cycle finding dfs, regarding only set edges
            memset(chk, -1, sizeof(long)*graph->num_vertices);
            bool cyclic = false;
            for (long u=0; u<graph->num_vertices; ++u)
            {
                if (chk[u] < 0)
                {
                    dfs_for_cyclic(u, -1, chk, original_list, set_edges, cyclic);
                    if (cyclic)
                        break;
                }
            }
            
            if (cyclic)
            {
                // TODO: remove me
                cout << "[chain_reduction] removing edge because it implies a cycle!" << endl;

                remove_edge(first_e);

                // TODO: remove me
                //cout << "\t\t\t\t now the graph has " << graph->num_edges << " edges:  ";
                //for (long e=0; e<graph->num_edges; ++e)
                //    cout << "#" << e << "(" << graph->s[e] << "," << graph->t[e] << ")  ";
                //cout << endl;
                //cout << "\t\t\t\t now the conflict graph has " << num_conflicts << " conflicts:   ";
                //for (long p=0; p<num_conflicts; ++p)
                //{
                //    pair<long,long> edges = conflicts.at(p);
                //    cout << "#" << edges.first << " vs #" << edges.second << "   ";
                //}
                //cout << endl;

                delete[] chk;
                delete[] lowest;
                delete[] inactive;
                delete[] set_edges;

                return 1;
            }
        }
        */
        
    }  // try new edge to start chain (could not derive infeasibility or unicity of solution)
    
    delete[] chk;
    delete[] lowest;
    delete[] inactive;
    delete[] set_edges;
    
    return 0;
}

// finds if "conflict closure" of fixing edge pairs would disconnect graph
long IO::preprocess_pairwise()
{
    long new_conflicts = 0;

    // 1. additional data structures for bridge finding dfs
    
    // original graph adjacency list
    vector< vector<long> > original_list;
    for (long i=0; i<graph->num_vertices; ++i)
        original_list.push_back(vector<long>());

    for (long e=0; e<graph->num_edges; ++e)
    {
        long i = graph->s[e];
        long j = graph->t[e];
    
        original_list[i].push_back(j);
        original_list[j].push_back(i);
    }
    
    // conflict graph adjacency list
    vector< vector<long> > conflict_list;

    for (long i=0; i<graph->num_edges; ++i)
        conflict_list.push_back(vector<long>());

    for (long p=0; p<num_conflicts; ++p)
    {
        pair<long,long> edges = conflicts.at(p);
        long e1 = edges.first;
        long e2 = edges.second;
        
        conflict_list[e1].push_back(e2);
        conflict_list[e2].push_back(e1);
    }
    
    // next are some auxiliary structures, which are overwritten several times

    // checks visited vertices (stores discovery time) in dfs
    long *chk = new long[graph->num_vertices];

    // earliest vertex seen in dfs
    long *lowest = new long[graph->num_vertices];
    
    // used for temporary selection or removal of edges
    bool *inactive = new bool[graph->num_edges];
    bool *set_edges = new bool[graph->num_edges];
    
    // 2. check connectivity fixing edges in chain, starting with first_e and second_e
    for (long first_e=0; first_e<graph->num_edges; ++first_e)
    {
        for (long second_e=first_e+1; second_e<graph->num_edges; ++second_e)
        {
            // two isolated vertices in the conflict graph cannot yield any reduction
            if (conflict_list[first_e].empty() && conflict_list[second_e].empty())
                continue;
            
            // conflict must not exist yet
            bool skip_pair = false;
            for (unsigned adj=0; adj<conflict_list[first_e].size(); ++adj)
                if (second_e == conflict_list[first_e].at(adj))
                {
                    skip_pair = true;
                    break;
                }
            if (skip_pair)
                continue;
        
            memset(set_edges, false, sizeof(bool)*graph->num_edges);
            memset(inactive, false, sizeof(bool)*graph->num_edges);
        
            // 3. temporarily select first_e and second_e, and reject conflicting ones
            long count_set_edges = 2;
            set_edges[first_e] = true;
            set_edges[second_e] = true;
            
            for (unsigned adj=0; adj<conflict_list[first_e].size(); ++adj)
                inactive[conflict_list[first_e].at(adj)] = true;
                
            for (unsigned adj=0; adj<conflict_list[second_e].size(); ++adj)
                inactive[conflict_list[second_e].at(adj)] = true;
        
            // 4. perform bridge finding dfs, disregarding inactive edges
            vector<long> bridges_idx;
            memset(chk, -1, sizeof(long)*graph->num_vertices);
            memset(lowest, -1, sizeof(long)*graph->num_vertices);
    
            long components = 0;
            long counter = -1;
        
            // shall execute dfs only once if G is connected
            for (long u=0; u<graph->num_vertices; ++u)
            {
                if (chk[u] < 0)
                {
                    ++components;
                    dfs_for_chain(u, counter, -1, bridges_idx, chk, lowest, original_list, set_edges, inactive);
                }
            }
    
            // 5. selection leads to a disconnected graph: include conflict (first_e, second_e), done!
            if (components > 1)
            {
                // TODO: remove me
                //cout << "\t\t graph would be disconnected! including conflict (" << first_e << "," << second_e << ")" << endl;
            
                add_conflict(first_e, second_e);
                ++new_conflicts;
            
                continue;
            }
            
            bool checking = true;

            // 6. selection imply bridges in the original graph, which shall be fixed
            while (bridges_idx.size() > 0 && checking)
            {
                // NOTE: would there be a preferable choice among bridges?
                long bridge = bridges_idx.back();

                // TODO: remove me
                if (inactive[bridge])
                    cout << endl << endl << "UNEXPECTED: selected bridge was inactive!" << endl << endl;
                if (set_edges[bridge])
                    cout << endl << endl << "UNEXPECTED: selected bridge was already set!" << endl << endl;
            
                // 7. temporary selection (conflicting edges could be inactive - ok!)
                ++count_set_edges;
                set_edges[bridge] = true;
                for (unsigned adj=0; adj<conflict_list[bridge].size(); ++adj)
                    inactive[conflict_list[bridge].at(adj)] = true;
            
                // 8. perform a new dfs, finding bridges and connectivity status
                bridges_idx.clear();
                memset(chk, -1, sizeof(long)*graph->num_vertices);
                memset(lowest, -1, sizeof(long)*graph->num_vertices);
            
                components = 0;
                counter = -1;
            
                for (long u=0; u<graph->num_vertices; ++u)
                {
                    if (chk[u] < 0)
                    {
                        ++components;
                        dfs_for_chain(u, counter, -1, bridges_idx, chk, lowest, original_list, set_edges, inactive);
                    }
                }
            
                // 9. selection leads to a disconnected graph: include conflict (first_e, second_e), done!
                if (components > 1)
                {
                    add_conflict(first_e, second_e);
                    ++new_conflicts;
                    
                    checking = false;
                }
                else
                {
                    // 10. is current graph a (spanning, as it is connected) tree? stop chain
                    if (count_set_edges == graph->num_vertices-1)
                    {
                        // stop chain as no more bridges will be checked
                        bridges_idx.clear();
                    
                        // check if conflict constraints are satisfied
                        bool feasible = true;
                        for (long p=0; p<num_conflicts; ++p)
                        {
                            pair<long,long> edges = conflicts.at(p);
                            long e1 = edges.first;
                            long e2 = edges.second;

                            // remove conflict if any of the edges was fixed
                            if (set_edges[e1] && set_edges[e2])
                            {
                                feasible = false;
                                break;
                            }
                        }
                    
                        // store if integer feasible solution to use as a primal bound
                        if (feasible)
                        {
                            // TO DO: determine objective and save if it's the new incumbent
                            // cout << "[pairwise_chain] found (primal) integer feasible solution " << endl;
                        }
                        else // no feasible spanning tree including this edge .: include conflict
                        {
                            // TO DO: remove me
                            //cout << "[pairwise_chain] including conflict because unique spanning tree is not feasible!" << endl;

                            add_conflict(first_e, second_e);
                            ++new_conflicts;
                            
                            checking = false;
                        }
                    }
                }
            
            } // fix another bridge (bridge_idx was updated properly)
            
            /*
            // if a cycle is implied, and first_e makes the problem infeasible .: just remove it
            if (count_set_edges >= 3)
            {
                // 11. perform cycle finding dfs, regarding only set edges
                memset(chk, -1, sizeof(long)*graph->num_vertices);
                bool cyclic = false;
                for (long u=0; u<graph->num_vertices; ++u)
                {
                    if (chk[u] < 0)
                    {
                        dfs_for_cyclic(u, -1, chk, original_list, set_edges, cyclic);
                        if (cyclic)
                            break;
                    }
                }
            
                if (cyclic)
                {
                    // TODO: remove me
                    cout << "[pairwise_chain] including conflict because it implies a cycle!!!" << endl;

                    add_conflict(first_e, second_e);
                    ++new_conflicts;
                    checking = false;
                }
            }
            */
        
        }  // second_e loop: try new pair to start chain (could not derive infeasibility or unicity of solution)
    
    } // first_e loop
    
    delete[] chk;
    delete[] lowest;
    delete[] inactive;
    delete[] set_edges;
    
    return new_conflicts;
}

// coordinates graph reduction and variable fixing tests
bool IO::preprocess()
{
    objective_offset = 0;
    long fixed_edges = 0;
    long new_conflicts = 0;
    bool updated;
    
    // TODO: remove me
    long chain_iterations = 0;
    
    do
    {
        updated = false;
        
        // 1. fix bridges and remove conflicting edges while possible
        long bridge_reduction = preprocess_bridges();
        fixed_edges += bridge_reduction;
        
        // graph became disconnected
        if (bridge_reduction < 0)
        {
            cout << endl << "infeasible instance: original graph does not have any spanning tree" << endl;
            this->problem_infeasible = true;
            return false;
        }

        // graph became a tree but has conflicting edges
        if (graph->num_edges == graph->num_vertices-1 && conflicts.size() > 0)
        {
            cout << "infeasible instance: original graph does not have any conflict-free spanning tree" << endl;
            this->problem_infeasible = true;
            return false;
        }

        // problem solved during preprocessing if there are no conflicts
        if (conflicts.empty())
        {
            cout << "preprocessing resulted in a trivial mst instance" << endl;
            this->problem_solved = true;
            return false;
        }
        
        // 2. try further reduction by verifying edge fixations in chain
        
        long chain_reduction = preprocess_chain();
        fixed_edges += chain_reduction;
        
        // TODO: remove me
        //cout << chain_reduction << " edge fixed with chain preprocessing" << endl;
        chain_iterations += chain_reduction;
        
        updated = (chain_reduction > 0);
        
        if (!updated)
        {
            // try generating new conflicts with pairwise chain processing
            long pairwise_chain = preprocess_pairwise();
            new_conflicts += pairwise_chain;
            
            if (pairwise_chain > 0)
                updated = true;
        }
    }
    while (updated);

    // update remaining data structures: conflict_graph_adj_list and graph->adj_list
    if (fixed_edges > 0)
    {
        graph->adj_list.clear();
        graph->adj_list.reserve(graph->num_vertices);
        graph->adj_list.insert(graph->adj_list.begin(), graph->num_vertices, list<long>() );

        // the edge list ("s" and "t" vectors) is updated correctly throughout the preprocessing algorithms
        for (long edge_idx=0; edge_idx<graph->num_edges; ++edge_idx)
        {
            long v1 = graph->s.at(edge_idx);
            long v2 = graph->t.at(edge_idx);
            
            graph->adj_list[v1].push_back(v2);
            graph->adj_list[v2].push_back(v1);
        }

        conflict_graph_adj_list.clear();
        conflict_graph_adj_list.reserve(graph->num_edges);
        conflict_graph_adj_list.insert(conflict_graph_adj_list.begin(), graph->num_edges, list<long>() );

        // the "conflicts" vector is correctly updated throughout the preprocessing algorithms
        for (vector< pair<long,long> >::iterator it = conflicts.begin(); it != conflicts.end(); ++it)
        {
            long e1_index = it->first;
            long e2_index = it->second;

            conflict_graph_adj_list[e1_index].push_back(e2_index);
            conflict_graph_adj_list[e2_index].push_back(e1_index);
        }
    }

    // TO DO: save the list of fixed_vars, and the map preprocessed_to_original_edge_idx
    
    cout << "preprocessing complete: "
         << fixed_edges << " fixed edges, "
         << new_conflicts << " new conflicts ";
    cout << "(new instance size: |V|=" << graph->num_vertices << ", |E|=" << graph->num_edges
         << ", |C|=" << num_conflicts << ")" << endl;

    // TO DO: remove me
    /*
    cout << endl << chain_iterations << " successfuly fixed edges by chain reduction" << endl;
    cout << endl << new_conflicts << " new conflicts introduced by pairwise chain" << endl << endl;
    
    cout << fixed_edges << " edges fixed by preprocessing algorithm"
         << " (objective offset = " << objective_offset << ")" << endl;
    
    cout << "reduced instance size: |V|=" << graph->num_vertices << ", |E|="
         << graph->num_edges << ", |C|=" << num_conflicts << endl;
    */
    
    return true;
}
