#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <map>
#include <set>
#include <sstream>
#include "graph.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "graph.h"
#include <time.h>
#include <queue>

using namespace std;
typedef double captype;

template <class InputIterator, class OutputIterator, class UnaryPredicate>
  void copy_if_v1 (InputIterator first, InputIterator last,
                          InputIterator in1, InputIterator in2,
                          OutputIterator result1, OutputIterator result2,
                          UnaryPredicate pred)
{
  while (first!=last) {
    if (pred(*first)) {
      *result1 = *in1;
      ++result1;
      *result2 = *in2;
      ++result2;
    }
    ++first;
    ++in1;
    ++in2;
  }
}

void read_nodes(const string& from, const string& to, vector<string>& nodes_fil, vector<string>& diffs_fil){

  ifstream ip("../_data/taxi-small/sunday-nodes.tsv");

  if(!ip.is_open()) std::cout << "ERROR: File Open" << '\n';
  string node;
  string diff;
  string group;
  string tmp;

  vector<string> nodes;
  vector<string> diffs;
  vector<string> groups;

  getline(ip, tmp);

  int i = 0;
  while(getline(ip, tmp)){
    i += 1;
    //if(i == 100) break;
    istringstream ss(tmp);
    getline(ss, node,'\t');
    getline(ss, diff,'\t');
    getline(ss, group,'\n');

    //cout << node << ' ' << diff << ' ' << group << ' ' << i << endl;
    nodes.push_back(node);
    diffs.push_back(diff);
    groups.push_back(group);
  }

  copy_if_v1(groups.begin(),
             groups.end(),
             nodes.begin(),
             diffs.begin(),
             std::inserter(nodes_fil, nodes_fil.begin()),
             std::inserter(diffs_fil, diffs_fil.begin()),
             [from,to](string& val){
  	return val >= from && val < to;
  });

  std::ofstream output_file("../example.txt");
  ostream_iterator<std::string> output_iterator(output_file, "\n");
  std::copy(nodes_fil.begin(), nodes_fil.end(), output_iterator);
  ip.close();
}

int read_edges(const vector<string>& nodes, vector<int>& srcs_fil_no_dup, vector<int>& trgs_fil_no_dup){

  ifstream ip("../_data/taxi-small/sunday-edges.tsv");

  if(!ip.is_open()) std::cout << "ERROR: File Open" << '\n';
  string src;
  string trg;
  string tmp;
  vector<string> srcs, trgs;

  getline(ip, tmp);

  int i = 0;
  while(getline(ip, tmp)){
    i += 1;
    //if(i == 100) break;
    istringstream ss(tmp);

    getline(ss, src, '\t');
    getline(ss, trg, '\n');

    srcs.push_back(src);
    trgs.push_back(trg);
    //std::cout << src << " ---- " << trg << '\n';
  }
  map<string, int> map_nodes;
  set<pair<int, int>> set_edges;
  for (i = 0; i < nodes.size(); i++){
    map_nodes[nodes[i]] = i;
  }
  vector<int> srcs_fil, trgs_fil;
  auto src_iter = srcs.begin();
  auto trg_iter = trgs.begin();
  while(src_iter != srcs.end() && trg_iter != trgs.end()){
    auto it1 = map_nodes.find(*src_iter);
    auto it2 = map_nodes.find(*trg_iter);

    if (it1 != map_nodes.end() && it2 != map_nodes.end()){
      if(it1->second <= it2->second){
        srcs_fil.push_back(it1->second);
        trgs_fil.push_back(it2->second);
      }
      else{
        srcs_fil.push_back(it2->second);
        trgs_fil.push_back(it1->second);
      }
    }
    src_iter++;
    trg_iter++;
  }


  for(i = 0; i < srcs_fil.size(); i++){
    pair<int, int> p(srcs_fil[i], trgs_fil[i]);
    if(set_edges.find(p) == set_edges.end()){
      set_edges.insert(p);
      srcs_fil_no_dup.push_back(srcs_fil[i]);
      trgs_fil_no_dup.push_back(trgs_fil[i]);
    }
  }

  ofstream out( "../test.txt" );
  for(i = 0; i < srcs_fil_no_dup.size(); i++){
    out << srcs_fil_no_dup[i] << ' ' << trgs_fil_no_dup[i] << endl;
    // cout << srcs_fil[i] << ' ' << trgs_fil[i] << endl;
    // cout << srcs_fil_no_dup[i] << ' ' << trgs_fil_no_dup[i] << endl;
  }
  out.close();

  ip.close();

  return set_edges.size();

}

/*-----------------------------------------------*/
int bfs(double** rGraph, int s, int t, int parent[], const int V) {
    // Create a visited array and mark all vertices as not visited
    bool *visited = new bool[V];
    memset(visited, 0, V*sizeof(visited[0]));

    // Create a queue, enqueue source vertex and mark source vertex
    // as visited
    queue <int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;

    // Standard BFS Loop
    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v = 0; v < V; v++) {
            if (visited[v] == false && rGraph[u][v] > 0) {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }

    // If we reached sink in BFS starting from source, then return
    // true, else false
    return (visited[t] == true);
}

// A DFS based function to find all reachable vertices from s.  The function
// marks visited[i] as true if i is reachable from s.  The initial values in
// visited[] must be false. We can also use BFS to find reachable vertices
void dfs(double** rGraph, int s, bool visited[], const int V) {
    visited[s] = true;
    for (int i = 0; i < V; i++)
       if (abs(rGraph[s][i]) > 1e-6 && !visited[i])
           dfs(rGraph, i, visited, V);
}

// Prints the minimum s-t cut
double** minCut(double** graph, int s, int t, bool *visited, const int V) {
    int u, v;
    double max_flow = 0;
    // Create a residual graph and fill the residual graph with
    // given capacities in the original graph as residual capacities
    // in residual graph
    double** rGraph = new double*[V];  // rGraph[i][j] indicates residual capacity of edge i-j
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

        // Find minimum residual capacity of the edges along the
        // path filled by BFS. Or we can say find the maximum flow
        // through the path found.
        double path_flow = INT_MAX;
        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            path_flow = min(path_flow, rGraph[u][v]);
        }

        // update residual capacities of the edges and reverse edges
        // along the path
        //printf("residual\n");
        for (v = t; v != s; v = parent[v]) {
            u = parent[v];
            //printf("%d -> %d\n", u, v);
            rGraph[u][v] -= path_flow;
            rGraph[v][u] += path_flow;
        }
        /*
        printf("residual graph printout (%d) \n", counter);
        for (u = 0; u < V; u++) {
            for (v = 0; v < V; v++) {
                 printf("%5.2f ", rGraph[u][v]);
            }
            printf("\n");
        }
        printf("\n");
        */
        counter++;
        max_flow += path_flow;
    }

    // Flow is maximum now, find vertices reachable from s
    memset(visited, false, V*sizeof(visited[0]));
    dfs(rGraph, s, visited, V);
    //for (int i = 0; i < V; i++) printf("visited: %d \n", visited[i]);

    // Print all edges that are from a reachable vertex to
    // non-reachable vertex in the original graph
    /*
    printf("Cut:\n");
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (visited[i] && !visited[j] && graph[i][j]) {
                cout << i << " -> " << j << endl;
            }
        }
    }
    */
    //printf("maxflow is %f \n", max_flow);
    return rGraph;
}
/*----------------------------------------------*/


void graph_tv(double *Y,// value of nodes
        int n, //number of nodes
        int m, // number of edges
        int* e1,
        int* e2,
//an array of edges of size m. There is an edge edges1[i] -> edges2[i]
        double lambda,
        float erreur)
{

    double lambda1 = lambda;
    int V = n+2;
    double *tr_cap = new double[n];
    unsigned int label[n];

    double **graph = new double*[V];
    for (int h = 0; h < V; h++){
      graph[h] = new double[V];
      for (int w = 0; w < V; w++)
        graph[h][w] = 0;
    }

    for (int i = 0; i < m; i++) {
      graph[e1[i]+2][e2[i]+2] = lambda1;
      graph[e2[i]+2][e1[i]+2] = lambda1;
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

    moy=1.25; //need to change later !!!
    moy=0.; num=0;
    for (int i = 0; i < n; i++) {
      moy += Y[i];
      alive[i]=1;
      num++;
    }
    moy /= (double) num;
    values[0]=moy;
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
        printf("iter %d\n", iter++);
        for (int u = 0; u < V; u++) {
            graph_old[u] = new double[V];
            for (int v = 0; v < V; v++) {
                 graph_old[u][v] = graph[u][v];
            }
        }
      /*printf("\n residual graph before min-cut: \n");
      for (int u = 0; u < V; u++) {
          for (int v = 0; v < V; v++) {
               printf("%5.2f ", graph[u][v]);
          }
          printf("\n");
      }
      printf("\n");
      */
      bool *visited = new bool[V];
      graph = minCut(graph, 0, 1, visited, V);

      /*
      printf("\n residual graph after min-cut: \n");
      for (int u = 0; u < V; u++) {
          for (int v = 0; v < V; v++) {
               printf("%f ", graph[u][v]);
          }
          printf("\n");
      }
      printf("\n");
      */

      memset(averages,0,nlab*sizeof(captype));
      memset(nums,0,nlab*sizeof(int));
      memset(nextlabel,0,nlab*sizeof(int));
      oldnlab=nlab;

        for (int i = 0; i < n; i++)
          if (alive[i]) {
            //printf("!!!!!!!!!!!!!! %d  \n",l);
            if (visited[i+2]==1) {
              l=nextlabel[label[i]];
              if (l==0) {
                l=(nextlabel[label[i]]=nlab);
                inactivelabel[l]=0;
                nlab++;
                averages[l]=0.; nums[l]=0;
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
              //printf("I am in ++ \n");
            } else { // what_segment(i) == sink
              l=label[i];
              averages[l] -= graph[i+2][1];
              nums[l]++;
              //	      i->tr_cap += delta;
              for (int j = 0; j < n; j++) // easy to be wrong !!!
                if (graph_old[i+2][j+2] != 0 && visited[j+2]) { graph[i+2][j+2]=0;}
              //printf("I am in -- \n");
            }
          }

          // tentative d'arret a precision zero
          // detection d'un label qui n'a pas ete coupe
          for (l=0;l<oldnlab;l++) if (!inactivelabel[l]) {
            if (nextlabel[l]==0) { averages[l]=0.; inactivelabel[l]=1;}
            else if (nums[l]==0) {
              inactivelabel[l]=inactivelabel[nextlabel[l]]=1;
              averages[nextlabel[l]]=0.;
            } else {
              averages[l] /= (double) nums[l];
              values[l] += averages[l];
            }
          } else averages[l]=0.;
          for (; l<nlab; l++) {
            averages[l] /= (double) nums[l];
            values[l]   += averages[l];
            //printf("!!!!!!!!!!!!!! %d, %d, %f \n",nums[l], l, averages[l]);
          }
        flagstop=0;
        //printf("\n residual graph after gtf!!!!!!!!!!!: \n");
        for (int i = 0; i < n; i++)
          if (alive[i]) {
            l = label[i];
            if (inactivelabel[l] ||
            (averages[l]<=erreur && averages[l]>=-erreur)) {
            if (visited[i+2]==1) graph[0][i+2] = 0; //!!!
            if (visited[i+2]==0) graph[i+2][1] = 0; //!!!
            alive[i] = 0; // noeud d�connect� � l'avenir
            inactivelabel[l]=1;
            } // end if
            else {
            flagstop=1; // on continue
            if (visited[i+2]==1) {
              graph[0][i+2] -= averages[l]; //!!!
              if(graph[0][i+2] < 0){
                double tmp = -graph[0][i+2];
                graph[0][i+2] = graph[i+2][1];
                graph[i+2][1] = tmp;
              }
            }
            if (visited[i+2]==0) {
              graph[i+2][1] += averages[l]; //!!!
              if(graph[i+2][1] < 0){
                double tmp = -graph[i+2][1];
                graph[i+2][1] = graph[0][i+2];
                graph[0][i+2] = tmp;
              }
            } // end else
          } // end for
        }
          /*
          printf("\n residual graph after gtf: \n");
          for (int u = 0; u < V; u++) {
              for (int v = 0; v < V; v++) {
                   printf("%5.2f ", graph[u][v]);
              }
              printf("\n");
          }
          printf("\n");
          */

        } while (flagstop);

        free(nextlabel);
        free(inactivelabel);
        free(nums);
        free(averages);

        for(int i = 0; i < n; i++)
          Y[i] = (double) values[label[i]];

}

void soft_thresh(double *Y, const double thresh, const int n){
    for(int i = 0; i < n; i++){
        double tmp = max(Y[i] - thresh, 0.0);
        Y[i] = tmp + min(Y[i]+thresh, 0.0);
    }
}


int main(){
    string start_time = "2011-06-26 12:00:00";
    string end_time = "2011-06-26 14:00:00";
    double lambda1 = 6;
    double lambda2 = 3;
    vector<string> nodes_fil, diffs_fil;
    read_nodes(start_time, end_time, nodes_fil, diffs_fil);

    int n = nodes_fil.size();
    vector<int> srcs_fil_no_dup, trgs_fil_no_dup;
    int m = read_edges(nodes_fil, srcs_fil_no_dup, trgs_fil_no_dup);

    double *Y = new double[n]; // nodes filled
    for(int i = 0; i < n; i++) Y[i] = stod(diffs_fil[i]);

    int *edges1, *edges2;
    edges1 = new int[m];
    edges2 = new int[m];

    for(int i = 0; i < m; i++){
      edges1[i] = srcs_fil_no_dup[i];
      edges2[i] = trgs_fil_no_dup[i];
    }

    cout << "Done! "<< "# of nodes: " << n << "; # of edges: " << m << endl;
    for(int i = 0; i < 30; i++){
      cout << Y[i] << ' ' << edges1[i] << ' ' << edges2[i] << endl;
    }

    clock_t t1, t2;
    t1 = clock();
	  graph_tv(Y, n, m, edges1, edges2, lambda1, 0.0);
    soft_thresh(Y, lambda2, n);
    t2 = clock();
    cout << "time is " << ((float)t2 - (float)t1) / CLOCKS_PER_SEC << endl;
    for(int i = 0;  i < 30; i++){
        printf("result is %f. \n", Y[i]);
    }
    ofstream out( "../output_min.txt" );
    //out.precision(3);
    for(int i = 0; i < n; i++){
        out << Y[i] << endl;
    }
    return 0;
}
