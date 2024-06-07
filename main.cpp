#include <iostream>
#include "src/graph.h"
#include "src/io.h"
#include "src/stack.h"
#include "src/performance.h"
#include "src/reductions.h"
#include "src/branchandreduce.h"
#include <signal.h>
#include <unistd.h>
#include <cstring>

/*
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
*/

// Global variables to store the best solution and its crossings
std::vector<Node*> bestSolution;
int bestCrossings;

// Signal handler for SIGTERM
void term(int signum) {
    // Output the best solution so far
    //std::cout << bestCrossings << std::endl;
    outputStandardOut(bestSolution);
    std::exit(signum); // Exit the program
}

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    // Register the signal handler for SIGTERM
    struct sigaction action;
    memset(&action, 0, sizeof(struct sigaction));
    action.sa_handler = term;
    sigaction(SIGTERM, &action, NULL);

    // Read graph
    Graph* g = readStandardIn();

    g->MedianHeuristic();
    bestSolution = g->getOrderNodes();
    bestCrossings = g->countCrossings();

    //std::cout << bestCrossings << ",";
    
    // Which reductions to use
    std::vector<general_reduction*> reductions;
    reductions.push_back(new ZeroEdge_reduction);
    reductions.push_back(new Complete_reduction);
    //reductions.push_back(new ZeroCrossings_reduction);
    reductions.push_back(new Twins_reduction);
    //reductions.push_back(new AlmostTwin_reduction);
    //reductions.push_back(new Domination_reduction);

    // Output solution
    int method2 = 2;
    bool fast = 1;

    TomAlvAlg(*g);
    int currentCrossings = g->countCrossings();
    if (currentCrossings < bestCrossings){
        bestSolution = g->getOrderNodes();
        bestCrossings = currentCrossings;
    }

    //std::cout << g->countCrossings() << ",";

    auto result = BranchAndReduce(g, reductions, 2, method2, fast);
    currentCrossings = g->countCrossings();
    if (currentCrossings < bestCrossings){
        bestSolution = result.first;
        bestCrossings = currentCrossings;
    }

    
    while(true){
        g->MedianHeuristic();
        auto result = BranchAndReduce(g, reductions, 3, method2, fast);
        int currentCrossings = g->countCrossings();
        if (currentCrossings < bestCrossings){
            bestSolution = result.first;
            bestCrossings = currentCrossings;
        }
    }
    
    
    //std::cout << bestCrossings << std::endl;
    //outputStandardOut(bestSolution);
}
