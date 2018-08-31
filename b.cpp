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


int main(){
    double lambda1 = 6;
    double lambda2 = 3;
    
    string e_file_name = "../_data/e.txt";
    string n_file_name = "../_data/n.txt";
    
    int *edges1, *edges2;
    double *Y;

    ifstream n_infile(n_file_name);
    ifstream e_infile(e_file_name);
    int tm1, tm2, j = 0, n = 0, m = 0;
    double tm3;  
    e_infile >> tm1 >> tm2;
    m = tm2;
    n = tm1;
    cout << m << n << "!!!!!!!!!!!!";
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
    ofstream out( "output_wang.txt" );
    //out.precision(3);
    for(int i = 0; i < n; i++){
        out << Y[i] << endl;
    }
    return 0;
}
