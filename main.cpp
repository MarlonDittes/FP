#include <iostream>
#include "graph.h"
#include "io.h"
#include "stack.h"
/*#include "performance.h"*/

constexpr int MODE_COUNT = 4;

int main(int argc, char* argv[]) {
    std::string graph_file = argv[1];
    std::cout << graph_file << std::endl;

    Graph* g = readGraph(graph_file);
    Graph verifier = *g;

    g->sortNeighbours();
    //g->printGraph();
    long crossing_count = g->countCrossings();
    std::cout << "number of crossings in default g: " << crossing_count << std::endl;

    long crossing_count_test = g->countCrossingsBranching();
    std::cout << "number of crossings in default g: " << crossing_count_test << std::endl;

    /* // Reduction Stuff
    g->Sorted_straight_line_reduction();
    g->printGraph();
    std::cout << "number of crossings in g: " << g->countCrossings() << std::endl;

    Twins twins = g->findTwins();
    std::cout << "twins size: " << twins.size() << std::endl;

    std::cout << "twins: " << std::endl;
    for (auto& twin: twins) {
        std::cout << twin.first->id << " - " << twin.second->id << std::endl;
    }

    std::cout << "crossings between 4, 5: " << g->countCrossingsForPair(1, 2) << std::endl;
    */

    // Testing Stuff
    // Brute Force (parallel)
    //auto resultBF = bruteForce(g);
    //auto resultBF = bruteForceParallel(*g);
    //std::cout << "exact crossings BRUTE FORCE: " << resultBF.second << std::endl;
    //outputOrder(resultBF.first, "../output.txt");



    // Greedy
    //auto resultGREEDY = g->Greedy();
    //std::cout << "crossings GREEDY: " << resultGREEDY.second << std::endl;

    // Median Heuristic
    //g->Median_Heuristic();
    //std::cout << "crossings MEDIAN_HEURISTIC: " << g->countCrossings() << std::endl;



    //g->printGraph();
    
    g->Partition();
    g->AP();
    g->printGraph();


    //Branch_and_Bound(g);
    

    if (g->verifier(verifier)) {
        std::cout << "Graph is valid"<<std::endl;
    }
    else {
        std::cout << "Graph is NOT valid" << std::endl;
    }

    // Testing on test set
    /*
    for (int mode = 1; mode<=MODE_COUNT; mode++){
        calculatePerformance("../tiny_test_set", "../performance_data/" + std::to_string(mode), mode);
        std::cout << "Mode " << mode << " done." << std::endl;
    } 
    */

    //readWeightedHyperGraph("../hypergraphs/4000_4_40_5.4.hgr", "../bipartitegraphs/4000_4_40_5.4.gr");

}

