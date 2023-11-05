#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <vector>
#include "graph.h"

Graph readGraphFile(std::string filename);
void printAdjacencyList(const std::vector<std::vector<int>>& adjacencyList);
bool verifier (Graph graph, std::vector<int> ordering);

#endif