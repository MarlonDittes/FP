#include <iostream>
#include "graph.h"
#include "io.h"
#include "stack.h"
#include "performance.h"
#include "reductions.h"

int main(int argc, char* argv[]) {
    // Read graph
    Graph* g = readStandardIn();

    //g->MedianHeuristic();
    //std::cout << g->countCrossingsMarlon() << ",";
    
    // Which reductions to use
    std::vector<general_reduction*> reductions;
    reductions.push_back(new ZeroEdge_reduction);
    reductions.push_back(new Complete_reduction);
    //reductions.push_back(new ZeroCrossings_reduction);
    reductions.push_back(new Twins_reduction);
    //reductions.push_back(new AlmostTwin_reduction);
    //reductions.push_back(new Domination_reduction);

    // Output solution
    int method1 = 1;
    int method2 = 2;
    bool fast = 1;
    auto result = BranchAndReduce(g, reductions, method1, method2, fast);

    std::cout << g->countCrossingsMarlon() << std::endl;
    //outputStandardOut(result.first);
}
