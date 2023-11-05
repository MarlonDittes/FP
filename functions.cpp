#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include "graph.h"
#include <set>

//AdjacencyList functionality
Graph readGraphFile(std::string filename){
    std::cout << "Reading Graph file..." << std::endl;
    std::ifstream inputFile(filename);
    std::vector<std::vector<int>> adjList;

    //Check if file can be opened
    if (!inputFile.is_open()){
        std::cout << "Failed to open " << filename << "." << std::endl;
        Graph graph;
        return graph;
    }

    std::string line;
    int n0, n1, m;
    int x,y;


    while (getline(inputFile, line)){
        //Skip empty lines or lines starting with 'c' (comments)
        if (line.empty() || line[0] == 'c'){
            continue;
        }

        std::istringstream iss(line);
        //p-line case
        if (line[0] == 'p'){
            std::string tmp;
            iss >> tmp >> tmp;
            iss >> n0 >> n1 >> m;
            
            adjList.resize(n0);
        } 
        //read edge case
        //we do -1 when reading so we can start at 0 on the vector
        else {
            iss >> x >> y;
            assert(1 <= x <= n0 && n0 + 1 <= y <= n0 + n1);
            x--;
            y--;

            adjList[x].push_back(y);
        }

    }
    Graph graph(adjList, n0, n1, m);
    return graph;
}

void printAdjacencyList(const std::vector<std::vector<int>>& adjacencyList) {
    for (int i = 0; i < adjacencyList.size(); ++i) {
        std::cout << "Vertex " << i << " is adjacent to: ";
        for (int j = 0; j < adjacencyList[i].size(); ++j) {
            std::cout << adjacencyList[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// OSCM tools
bool verifier (Graph graph, std::vector<int> ordering){
    std::set<int> uniqueElems;  // To store unique numbers, no node allowed to occur twice

    for (int elem : ordering){
        //Assume we still need to go back to +1 notation
        elem++;

        //Case: node is not element of B (free order)
        if (elem <= graph.n0 || elem > graph.n0 + graph.n1){
            return false;
        }

        //Check if node already occured
        if (uniqueElems.find(elem) != uniqueElems.end()) {
            return false;
        }

        //Add number to set
        uniqueElems.insert(elem);
    }

    //All elems within range and unique
    return true;
}

std::vector<int> bruteForce (Graph graph){

}

