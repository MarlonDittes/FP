#ifndef BRANCHANDREDUCE_H
#define BRANCHANDREDUCE_H

#include "graph.h"
#include "../src_tomalv/heuristic_algorithm.h"
#include "../src_tomalv/heuristic_graph.h"
#include "../src_tomalv/median_algorithm.h"

#include <cstdlib> // For rand() and srand()
#include <ctime>   // For time()
#include <memory>

int calculateSpan(Node* node);
std::unique_ptr<Graph> createGraphByPartition(std::unique_ptr<Graph> g, std::vector<Node*> partition);
void TomAlvAlg(Graph& g);
std::pair<std::vector<Node*>, long> branching(Graph* g, std::vector<std::unique_ptr<general_reduction>>& reductionTypes, int method1, int method2, bool fast);
std::pair<std::vector<Node*>, long> BranchAndReduce(Graph* g, std::vector<std::unique_ptr<general_reduction>>& reductionTypes, int method1, int method2, bool fast);
std::pair<std::vector<Node*>, long>  ExactSolution(Graph& g);

#endif //BRANCHANDREDUCE_H