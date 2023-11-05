#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>
#include <vector>

class Graph{
    public:
        std::vector<std::vector<int>> adjList;
        std::vector<int> fixedOrder;
        //std::vector<int> freeOrder;

        int n0, n1, m;

        Graph():
            n0(-1){};
        Graph(std::vector<std::vector<int>> AdjList, int N0, int N1, int M):
            adjList(AdjList), n0(N0), n1(N1), m(M){};
};

#endif