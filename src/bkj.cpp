#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <map>
#include <set>
#include <sstream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <queue>

#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace std;
typedef double captype;

// comment out to use FFA
#define push_relabel


/*-----------------------------------------------*/
int bfs(double** rGraph, int s, int t, int parent[], const int V) {
    // Create a visited array and mark all vertices as not visited
    bool *visited = new bool[V];
    memset(visited, 0, V*sizeof(visited[0]));

    queue <int> q;
    q.push(s);
    visited[s] = true;
    parent[s]  = -1;

    // Standard BFS Loop
    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v = V - 1; v >= 0; v--) {
            if (visited[v] == false && rGraph[u][v] > 0) {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }

    return (visited[t] == true);
}

// A DFS based function to find all reachable vertices from s.  The function
// marks visited[i] as true if i is reachable from s.  The initial values in
// visited[] must be false. We can also use BFS to find reachable vertices
void dfs(double** rGraph, int s, bool visited[], const int V) {
    visited[s] = true;
    for (int i = V - 1; i >= 0; i--) {
       if (rGraph[s][i] != 0 && !visited[i])
           dfs(rGraph, i, visited, V);
    }
}

// Prints the minimum s-t cut
void ffa_mincut(double** graph, int s, int t, bool *visited, const int V) {
    int u, v;
    double max_flow = 0;
    double** rGraph = new double*[V];
    for (u = 0; u < V; u++) {
        rGraph[u] = new double[V];
        for (v = 0; v < V; v++) {
             rGraph[u][v] = graph[u][v];
        }
    }

    int *parent = new int[V];  // This array is filled by BFS and to store path

    // Augment the flow while there is a path from source to sink
    int counter = 0;
    while (bfs(rGraph, s, t, parent, V)) {
        double path_flow = INT_MAX;
        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            path_flow = min(path_flow, rGraph[u][v]);
        }
        // std::cerr << "path_flow=" << path_flow << std::endl;

        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }
        counter++;
        max_flow += path_flow;
    }

    memset(visited, false, V*sizeof(visited[0]));
    dfs(rGraph, s, visited, V);

    double z = 0;
    for(int i = 0; i < V; i++) {
        for(int j = 0; j < V; j++) {
            // std::cerr << graph[i][j] << ":" << rGraph[i][j] << " ";
            graph[i][j] = rGraph[i][j];
            z += rGraph[i][j];
        }
        // std::cerr << std::endl;
    }
    std::cerr << "ffa_max_flow=" << max_flow << " | z=" << z << std::endl;
}

void pr_mincut(double** graph, int s, int t, bool *visited, const int naug_nodes) {

    using namespace boost;

    typedef adjacency_list_traits < vecS, vecS, directedS > Traits;
    typedef adjacency_list < vecS, vecS, directedS,
        property < vertex_name_t, std::string >,
        property < edge_capacity_t, double,
        property < edge_residual_capacity_t, double,
        property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;

    Graph boost_graph;

    typename property_map < Graph, edge_capacity_t >::type
        capacity = get(edge_capacity, boost_graph);

    typename property_map < Graph, edge_reverse_t >::type
        rev = get(edge_reverse, boost_graph);

    typename property_map < Graph, edge_residual_capacity_t >::type
        residual_capacity = get(edge_residual_capacity, boost_graph);

    std::vector<Traits::vertex_descriptor> verts;
    for (int v = 0; v < naug_nodes; ++v) {
        verts.push_back(add_vertex(boost_graph));
    }

    Traits::vertex_descriptor source = verts[s];
    Traits::vertex_descriptor sink   = verts[t];

    for(int i = 0; i < naug_nodes; i++) {
        for(int j = 0; j < naug_nodes; j++) {
            if(graph[i][j] != 0) {
                Traits::edge_descriptor e1, e2;
                bool in1, in2;
                tie(e1, in1) = add_edge(verts[i], verts[j], boost_graph);
                tie(e2, in2) = add_edge(verts[j], verts[i], boost_graph);
                if (!in1 || !in2){
                    std::cerr << "error" << std::endl;
                }
                capacity[e1] = graph[i][j];
                capacity[e2] = 0;
                rev[e1] = e2;
                rev[e2] = e1;
            }
        }
    }

    double maxflow = push_relabel_max_flow(boost_graph, source, sink);

    std::vector<std::vector<double>> boost_resid;
    boost_resid.resize(naug_nodes);
    for (int i = 0; i < naug_nodes; i++) {
        boost_resid[i].resize(naug_nodes, 0.0);
        for(int j = 0; j < naug_nodes; j++) {
            boost_resid[i][j] = 0.0;
        }
    }

    typename graph_traits<Graph>::vertex_iterator u_it, u_end;
    typename graph_traits<Graph>::out_edge_iterator e_it, e_end;
    for (tie(u_it, u_end) = vertices(boost_graph); u_it != u_end; ++u_it){
        for (tie(e_it, e_end) = out_edges(*u_it, boost_graph); e_it != e_end; ++e_it){
            if (capacity[*e_it] > 0){
                int t = target(*e_it, boost_graph);
                boost_resid[*u_it][t] = residual_capacity[*e_it];
                boost_resid[t][*u_it] = (capacity[*e_it] - residual_capacity[*e_it]); // Added this
            }
        }
    }

    double z = 0;
    for(int i = 0; i < naug_nodes; i++) {
        for(int j = 0; j < naug_nodes; j++) {
            graph[i][j] = boost_resid[i][j];
            z += boost_resid[i][j];
        }
    }

    memset(visited, false, naug_nodes * sizeof(visited[0]));
    dfs(graph, s, visited, naug_nodes);

    std::cerr << "pr_max_flow=" << maxflow << " | z=" << z << std::endl;
}



void graph_tv(double *Y,// value of nodes
        int n, //number of nodes
        int m, // number of edges
        int* e1,
        int* e2,
//an array of edges of size m. There is an edge edges1[i] -> edges2[i]
        double lambda,
        float erreur){

    double lambda1 = lambda;
    int V = n+2;
    double *tr_cap = new double[n];
    unsigned int label[n];

    double **graph = new double*[V];
    for (int h = 0; h < V; h++){
        graph[h] = new double[V];
        for (int w = 0; w < V; w++) graph[h][w] = 0;
    }

    // This destroys all direction information?
    for (int edge_idx = 0; edge_idx < m; edge_idx++) {
        graph[e1[edge_idx]+2][e2[edge_idx]+2] = lambda1;
        graph[e2[edge_idx]+2][e1[edge_idx]+2] = lambda1;
    }

    // normalization as did in parametric maxflow
    #define Alloc_Size 1024
    unsigned int l,*nextlabel,nlab,oldnlab,*nums;
    unsigned char flagstop, *inactivelabel, *alive;
    double *averages;
    captype *values;
    unsigned int maxlabel; // the size of the array "values

    maxlabel=Alloc_Size;
    nextlabel= (unsigned int *) malloc(sizeof(unsigned int)*maxlabel);
    inactivelabel= (unsigned char *) malloc(sizeof(unsigned char)*maxlabel);
    nums = (unsigned int *) malloc(sizeof(unsigned int)*maxlabel);
    averages = (captype *) malloc(sizeof(captype)*maxlabel);
    values = (captype *) malloc(sizeof(captype)*maxlabel);
    alive = (unsigned char *) malloc(sizeof(unsigned char)*n);

    double moy;
    int num;

    nlab=1;
    nextlabel[0]=0;
    memset(inactivelabel,0,sizeof(unsigned char)*maxlabel);

    // moy=1.25; //need to change later !!!
    moy=0.; num=0;
    for (int i = 0; i < n; i++) {
        moy += Y[i];
        alive[i]=1;
        num++;
    }

    moy /= (double) num;
    values[0] = moy;

    for (int i = 0; i < n; i++) {
        double tem = Y[i] - moy;
        if(tem > 0) graph[0][i+2] = tem;
        if(tem < 0) graph[i+2][1] = -tem;
        label[i]=0;
    }
    int iter = 0;

    /* ----------------------------------------- */
    double** graph_old = new double*[V];  // rGraph[i][j] indicates residual capacity of edge i-j

    do {
        for (int u = 0; u < V; u++) {
            graph_old[u] = new double[V];
            for (int v = 0; v < V; v++) {
                graph_old[u][v] = graph[u][v];
            }
        }

        bool *visited = new bool[V];

#ifdef push_relabel
        pr_mincut(graph, 0, 1, visited, V);
#else
        ffa_mincut(graph, 0, 1, visited, V);
#endif

        memset(averages,  0, nlab * sizeof(captype));
        memset(nums,      0, nlab * sizeof(int));
        memset(nextlabel, 0, nlab * sizeof(int));
        oldnlab = nlab;

        for (int i = 0; i < n; i++) {
            if (alive[i]) {
                if (visited[i+2] == 1) {
                    l = nextlabel[label[i]];
                    if (l==0) {
                        l=(nextlabel[label[i]]=nlab);
                        inactivelabel[l]=0;
                        nlab++;
                        averages[l]=0.;
                        nums[l]=0;
                        nextlabel[l]=0;
                        values[l]=values[label[i]];
                        if (nlab==maxlabel) {
                            maxlabel+=Alloc_Size;
                            inactivelabel= (unsigned char *) realloc(inactivelabel,sizeof(unsigned char)*maxlabel);
                            nextlabel= (unsigned int *) realloc(nextlabel,sizeof(unsigned int)*maxlabel);
                            nums = (unsigned int *) realloc(nums,sizeof(unsigned int)*maxlabel);
                            averages = (captype *) realloc(averages,sizeof(captype)*maxlabel);
                            values = (captype *) realloc(values,sizeof(captype)*maxlabel);
                        }
                    } // end l == 0
                    label[i] = l;
                    averages[l] += graph[0][i+2]; // might be wrong!!!
                    nums[l]++;
                } else { // what_segment(i) == sink
                    l = label[i];
                    averages[l] -= graph[i+2][1];
                    nums[l]++;
                  for (int j = 0; j < n; j++) // easy to be wrong !!!
                    if (graph_old[i+2][j+2] != 0 && visited[j+2]) { graph[i+2][j+2]=0;}
                }
            }
        }

        for (l=0;l<oldnlab;l++) {
            if (!inactivelabel[l]) {
                if (nextlabel[l]==0) {
                    averages[l]=0.; inactivelabel[l]=1;
                } else if (nums[l]==0) {
                    inactivelabel[l]=inactivelabel[nextlabel[l]]=1;
                    averages[nextlabel[l]]=0.;
                } else {
                    averages[l] /= (double) nums[l];
                    values[l] += averages[l];
                }
            } else {
                averages[l]=0.;
            }
        }
        for (; l<nlab; l++) {
            averages[l] /= (double) nums[l];
            values[l]   += averages[l];
        }

        flagstop=0;

        for (int i = 0; i < n; i++) {
            if (alive[i]) {
                l = label[i];
                if (inactivelabel[l] || (averages[l]<=erreur && averages[l]>=-erreur)) {
                      if (visited[i+2]==1) graph[0][i+2] = 0; //!!!
                      if (visited[i+2]==0) graph[i+2][1] = 0; //!!!
                      alive[i] = 0;
                      inactivelabel[l]=1;
                } else {
                    flagstop=1; // on continue
                    if (visited[i+2] == 1) {
                        graph[0][i+2] -= averages[l]; //!!!
                        if(graph[0][i+2] < 0){
                            double tmp = -graph[0][i+2];
                            graph[0][i+2] = graph[i+2][1];
                            graph[i+2][1] = tmp;
                        }
                    }
                    if (visited[i+2] == 0) {
                        graph[i+2][1] += averages[l]; //!!!
                        if(graph[i+2][1] < 0){
                            double tmp = -graph[i+2][1];
                            graph[i+2][1] = graph[0][i+2];
                            graph[0][i+2] = tmp;
                        }
                    }
                }
            }
        }
    } while (flagstop);

   free(nextlabel);
   free(inactivelabel);
   free(nums);
   free(averages);
   for(int i = 0; i < n; i++) Y[i] = (double) values[label[i]];
}

void soft_thresh(double *Y, const double thresh, const int n){
    for(int i = 0; i < n; i++){
        Y[i] = max(Y[i] - thresh, 0.0) + min(Y[i] + thresh, 0.0);
    }
}

int main(int argc, char** argv){
    double lambda1 = 6;
    double lambda2 = 6;

    string n_file_name = argv[1];
    ifstream n_infile(n_file_name);

    string e_file_name = argv[2];
    ifstream e_infile(e_file_name);

    int num_nodes;
    int num_edges;
    e_infile >> num_nodes >> num_edges;

    // Read edge data
    int idx = 0;
    int src, dst;
    int* srcs = new int[num_edges];
    int* dsts = new int[num_edges];
    while(e_infile >> src >> dst){
      srcs[idx] = src;
      dsts[idx] = dst;
      idx++;
    }

    // Read node data
    idx = 0;
    double val;
    double* vals = new double[num_nodes];
    while(n_infile >> val){
      vals[idx] = val;
      idx++;
    }

    graph_tv(vals, num_nodes, num_edges, srcs, dsts, lambda1, 0.0);
    soft_thresh(vals, lambda2, num_nodes);

    for(int i = 0; i < num_nodes; i++) {
        cout << vals[i] << endl;
    }

    return 0;
}