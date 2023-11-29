#include <iostream>
#include "graph.h"
#include "io.h"
#include "stack.h"

int main(int argc, char* argv[]) {
    std::string graph_file = argv[1];
    std::cout << graph_file << std::endl;

    Graph* g = readGraph(graph_file);
    Graph verifier = *g;

    g->sortMovableNodes();
    std::cout << "sorted" << std::endl;

    g->printGraph();
    std::cout << "count crossings" << std::endl;
    std::cout << "number of crossings in g: " << g->countCrossings() << std::endl;

   /* Twins twins = g->findTwins();
    std::cout << "twins size: " << twins.size() << std::endl;

    std::cout << "twins: " << std::endl;
    for (auto& twin: twins) {
        std::cout << twin.first->id << " - " << twin.second->id << std::endl;
    }

    std::cout << "crossings between 4, 5: " << g->countCrossingsForPair(1, 2) << std::endl;
*/
    /*auto result = bruteForce(*g);
    std::cout << "crossings BRUTE FORCE: " << result.second << std::endl;
    //std::cout << "crossings GREEDY: " << g->Greedy().second << std::endl;

    //g->Median_Heuristic();
    //std::cout << "number of crossings in g: " << g->countCrossings() << std::endl;

    g->Partition();
    g->printGraph();
    g->AP();


    if (g->verifier(verifier)) {
        std::cout << "Graph is valid"<<std::endl;
    }
    else {
        std::cout << "Graph is NOT valid" << std::endl;
    }

    outputOrder(result.first, "../output.txt");
*/
    

readHyperGraph("../hypergraph.txt", "../graph2.txt");
}

