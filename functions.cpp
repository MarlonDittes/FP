#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>

//AdjacencyList functionality
std::vector<std::vector<int>> readGraphFile(std::string filename){
    std::cout << "Reading Graph file..." << std::endl;
    std::ifstream inputFile(filename);
    std::vector<std::vector<int>> adjList;

    //Check if file can be opened
    if (!inputFile.is_open()){
        std::cout << "Failed to open " << filename << "." << std::endl;
        return adjList;
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

    return adjList;
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

void verifier (const std::vector<std::vector<int>>& adjacencyList, std::vector<int> ordering){
    
}