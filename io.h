//INPUT OUTPUT FUNCTIONS

#ifndef IO_H
#define IO_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "graph.h"

Graph* readGraph(std::string graph_file);
void outputOrder(std::vector<Node*> order, std::string output);

void writeHyperGraphToBipartiteGraph(int n0, int n1, int m, std::vector<std::pair<int, int>> edges, std::string output);
void readHyperGraph(std::string hypergraph_file, std::string output);
void readWeightedHyperGraph(std::string hypergraph_file, std::string output);

#endif //IO_H