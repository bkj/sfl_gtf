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

using namespace std;

void graph_tv(double *Y,// value of nodes
        int n, //number of nodes
        int m, // number of edges
        int* edges1,
        int* edges2,
//an array of edges of size m. There is an edge edges1[i] -> edges2[i]
        double lambda,
        float erreur)
{
        int i;
    double dval,dr,drd;

    dr=lambda;

    Graph::node_id * nodes, no;

    if (n <=2 || m <= 2)
        { fprintf(stderr,"error: bad size\n"); exit(0); }

    nodes = (Graph::node_id *) malloc(n*sizeof(Graph::node_id));
    Graph *BKG = new Graph();
    for (i=0;i<n;i++) {
        no=nodes[i]=BKG->add_node();
        BKG->set_tweights(no,Y[i],0.);
    }
    for (i=0;i<m;i++) {
        BKG->add_edge(nodes[edges1[i]],nodes[edges2[i]],dr,dr);
    }
    BKG->dyadicparametricTV(erreur);
    for (i=0;i<n;i++) Y[i]=BKG->what_value(nodes[i]);
         delete BKG; free(nodes);
}

void soft_thresh(double *Y, const double thresh, const int n){
    for(int i = 0; i < n; i++){
        double tmp = max(Y[i] - thresh, 0.0);
        Y[i] = tmp + min(Y[i]+thresh, 0.0);
    }
}


int main(int argc, char** argv){
  double lambda1 = 6;
  double lambda2 = 3;
  int *edges1, *edges2;
  bool big_graph = false;
  int m, n;
  double *Y;

  string e_file_name = "./_data/e";
  string n_file_name = "./_data/n";

  ifstream n_infile(n_file_name);
  ifstream e_infile(e_file_name);
  int tm1, tm2, j = 0;
  double tm3;
  e_infile >> tm1 >> tm2;
  m = tm2;
  n = tm1;
  Y = new double[n]; // nodes filled
  edges1 = new int[m];
  edges2 = new int[m];
  while(e_infile >> tm1 >> tm2){
    edges1[j] = tm1;
    edges2[j++] = tm2;
  }

  j = 0;
  while(n_infile >> tm3){
    Y[j++] = tm3;
  }

  //cout << "Done! "<< "# of nodes: " << n << "; # of edges: " << m << endl;

  clock_t t1, t2;
  t1 = clock();
  graph_tv(Y, n, m, edges1, edges2, lambda1, 0.0);
  soft_thresh(Y, lambda2, n);
  t2 = clock();
  //cout << "time is " << ((float)t2 - (float)t1) / CLOCKS_PER_SEC << endl;
  for(int i = 0; i < n; i++){
    cout << Y[i] << endl;
  }
  return 0;
}
