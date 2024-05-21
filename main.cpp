#include <iostream>
#include "graph.h"
#include "io.h"
#include "stack.h"
#include "performance.h"
#include "reductions.h"

#include "src_henning/src/definitions.h"
#include "src_henning/src/macros.h"
#include "src_henning/src/misc.h"
#include "src_henning/src/graph_hen.h"
#include "src_henning/src/solver_bf.h"
#include "src_henning/src/exhaustive_solver.h"
#include "src_henning/src/partitioner.h"
#include "src_henning/src/solver.h"
#include "src_henning/src/useless_reducer.h"
#include "src_henning/src/front_back_reducer.h"

int main(int argc, char* argv[]) {
     // Check if there are enough command-line arguments
    //if (argc < 6) {
    //    std::cerr << "Usage: " << argv[0] << " <graph_file> <method1> <method2> <fast> <almost_twin>" << std::endl;
    //    std::cerr << "Types: string int int bool bool" << std::endl;
    //    return 1;
    //}

    std::string graph_file = argv[1];
    
    /*int method1 = std::stoi(argv[2]);
    int method2 = std::stoi(argv[3]);
    
    int boolArgInt = std::stoi(argv[4]);
    bool fast = (boolArgInt != 0);

    boolArgInt = std::stoi(argv[5]);
    bool almost = (boolArgInt != 0);*/

    //Random graph
    /*
    int sizeX = 6;
    int sizeY = 6;
    std::string output = "../bipartitegraphs/faulty.gr";
    auto edges = generateBipartiteGraph(sizeX, sizeY);
    writeGraphToBipartiteGraph(sizeX, sizeY, edges.size(), edges, output);
    */

    Graph* g = readGraph(graph_file);
    Graph verifier = *g;

    ExactSolution(g);

    //g->printGraph();
    //std::cout << std::endl;

    long crossing_count = g->countCrossingsMarlon();
    std::cout << "number of crossings in default g: " << crossing_count << std::endl;

    /*
    Graph* copy = readGraph(graph_file);
    copy->BarycenterHeuristicMarlon();
    crossing_count = copy->countCrossingsMarlon();
    std::cout << "number of crossings with BarycenterHeuristic g: " << crossing_count << std::endl;
    */

    //Graph* copy2 = readGraph(graph_file);
    //copy2->MedianHeuristic();
    //crossing_count = copy2->countCrossingsMarlon();
    //std::cout << "number of crossings with MedianHeuristic g: " << crossing_count << std::endl;

    //std::cout << std::endl;
    //
    //// Which reductions to use
    //std::vector<general_reduction*> reductions;
    //reductions.push_back(new ZeroEdge_reduction);
    //reductions.push_back(new Complete_reduction);
    ////reductions.push_back(new ZeroCrossings_reduction);
    //reductions.push_back(new Twins_reduction);
    //if (almost){
    //    reductions.push_back(new AlmostTwin_reduction);
    //}
    ////reductions.push_back(new Domination_reduction);

    //auto result = BranchAndReduce(g, reductions, method1, method2, fast);
    ////outputOrder(result.first, "../output2.txt");

    //if (!fast){
    //    if (result.second != g->countCrossingsMarlon()){
    //    std::cout << "Crossings didn't sum up!" << std::endl;
    //    } else {
    //        std::cout << "Crossings BranchAndReduce: " << result.second << std::endl;
    //    }
    //} else{
    //    std::cout << "Crossings BranchAndReduce: " << g->countCrossingsMarlon() << std::endl;
    //}

    //if (g->verifier(verifier)) {
    //    std::cout << "Graph is valid" << std::endl;
    //}
    //else {
    //    std::cout << "Graph is NOT valid" << std::endl;
    //}

    //for (auto& reduct : reductions){
    //    std::cout << reduct->get_reduction_type() << " was used " << reduct->usage_count << " times." << std::endl;
    //}   

    //Testing on test set

    /*
    std::vector<int> modes;
    //modes.push_back(3);
    modes.push_back(4);


    for (auto mode : modes){
        calculatePerformance("../exact-public", "../performance_data/" + std::to_string(mode), mode, method1, method2, fast, almost);
        std::cout << "Mode " << mode << " done." << std::endl;
    }
    */
    
}
