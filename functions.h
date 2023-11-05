#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <iostream>
#include <vector>

std::vector<std::vector<int>> readGraphFile(std::string filename);
void printAdjacencyList(const std::vector<std::vector<int>>& adjacencyList);

#endif