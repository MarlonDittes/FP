#include <iostream>
#include "graph.h"
#include "io.h"
#include "stack.h"
#include "performance.h"
#include "reductions.h"

int main(int argc, char* argv[]) {
    // Read graph
    Graph* g = readStandardIn();

    std::vector<Node*> bestSolution;
    int bestCrossings;

    g->MedianHeuristic();

    bestSolution = g->getOrderNodes();
    bestCrossings = g->countCrossingsMarlon();

    std::cout << bestCrossings << ",";
    
    // Which reductions to use
    std::vector<general_reduction*> reductions;
    reductions.push_back(new ZeroEdge_reduction);
    reductions.push_back(new Complete_reduction);
    //reductions.push_back(new ZeroCrossings_reduction);
    reductions.push_back(new Twins_reduction);
    //reductions.push_back(new AlmostTwin_reduction);
    //reductions.push_back(new Domination_reduction);

    // Output solution
    int method1 = 0;
    int method2 = 2;
    bool fast = 1;

    std::vector<int> method1_options = {0,1,2,3};
    for (auto& option : method1_options){
        auto result = BranchAndReduce(g, reductions, option, method2, fast);
        int currentCrossings = g->countCrossingsMarlon();
        std::cout << currentCrossings << ",";
        if (currentCrossings < bestCrossings){
            bestSolution = result.first;
            bestCrossings = currentCrossings;
        }
    }
    std::cout << std::endl;
    //outputStandardOut(result.first);
}
