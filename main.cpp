#include <iostream>
#include "graph.h"
#include "io.h"
#include "stack.h"
#include "performance.h"
#include "reductions.h"
#include <signal.h>
#include <unistd.h>
#include <cstring>

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
    bestCrossings = g->countCrossingsMarlon();

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

    std::vector<int> method1_options = {0,1,2};
    for (auto& option : method1_options){
        auto result = BranchAndReduce(g, reductions, option, method2, fast);
        int currentCrossings = g->countCrossingsMarlon();
        if (currentCrossings < bestCrossings){
            bestSolution = result.first;
            bestCrossings = currentCrossings;
        }
    }

    //std::cout << bestCrossings << std::endl;
    outputStandardOut(bestSolution);
}
