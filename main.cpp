#include <iostream>
#include "graph.h"
#include "io.h"
#include "stack.h"
#include "performance.h"
#include "reductions.h"

int main(int argc, char* argv[]) {
    std::string graph_file = argv[1];

    //Random graph
    /*
    int sizeX = 100;
    int sizeY = 100;
    std::string output = "../.Run/graph.txt";
    auto edges = generateBipartiteGraph(sizeX, sizeY);
    writeGraphToBipartiteGraph(sizeX, sizeY, edges.size(), edges, output);
    */

    Graph* g = readGraph(graph_file);
    Graph verifier = *g;

    //g->printGraph();
    //std::cout << std::endl;

    long crossing_count = g->countCrossingsMarlon();
    std::cout << "number of crossings in default g: " << crossing_count << std::endl;

    Graph* copy = readGraph(graph_file);
    copy->BarycenterHeuristicMarlon();
    crossing_count = copy->countCrossingsMarlon();
    std::cout << "number of crossings with BarycenterHeuristic g: " << crossing_count << std::endl;

    Graph* copy2 = readGraph(graph_file);
    copy2->MedianHeuristic();
    crossing_count = copy2->countCrossingsMarlon();
    std::cout << "number of crossings with MedianHeuristic g: " << crossing_count << std::endl;

    std::cout << std::endl;
    
    // Which reductions to use
    std::vector<general_reduction*> reductions;
    reductions.push_back(new ZeroEdge_reduction);
    reductions.push_back(new Complete_reduction);
    reductions.push_back(new ZeroCrossings_reduction);
    reductions.push_back(new Twins_reduction);

    auto result = BranchAndReduce(g, reductions);
    //outputOrder(result.first, "../output.txt");

    if (result.second != g->countCrossingsMarlon()){
        std::cout << "Crossings didn't sum up!" << std::endl;
    } else {
        std::cout << "Crossings BranchAndReduce: " << result.second << std::endl;
    }

    if (g->verifier(verifier)) {
        std::cout << "Graph is valid" << std::endl;
    }
    else {
        std::cout << "Graph is NOT valid" << std::endl;
    }

    //Testing on test set

    /*
    std::vector<int> modes;
    modes.push_back(1);
    modes.push_back(2);
    modes.push_back(3);
    modes.push_back(4);
    modes.push_back(5);


    for (auto mode : modes){
        calculatePerformance("../tiny_test_set", "../performance_data/" + std::to_string(mode), mode);
        std::cout << "Mode " << mode << " done." << std::endl;
    }
    */
}
