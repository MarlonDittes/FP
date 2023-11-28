//INPUT OUTPUT FUNCTIONS

#ifndef IO_H
#define IO_H
#include <iostream>
#include <fstream>
#include <sstream>
#include "graph.h"

Graph* readGraph(std::string graph_file);
void outputOrder(std::vector<Node*> order, std::string output);

#endif //IO_H