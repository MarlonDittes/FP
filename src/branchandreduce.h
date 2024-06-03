#ifndef BRANCHANDREDUCE_H
#define BRANCHANDREDUCE_H

#include "graph.h"
#include "../TomAlv/heuristic_algorithm.h"
#include "../TomAlv/heuristic_graph.h"
#include "../TomAlv/median_algorithm.h"

#include <cstdlib> // For rand() and srand()
#include <ctime>   // For time()

int calculateSpan(Node* node);
Graph* createGraphByPartition(Graph* g, std::vector<Node*> partition);
void TomAlvAlg(Graph& g);
std::pair<std::vector<Node*>, long> branching(Graph* g, std::vector<general_reduction*> reductionTypes, int method1, int method2, bool fast);
std::pair<std::vector<Node*>, long> BranchAndReduce(Graph* g, std::vector<general_reduction*> reductionTypes, int method1, int method2, bool fast);
std::pair<std::vector<Node*>, long>  ExactSolution(Graph& g);

#endif //BRANCHANDREDUCE_H