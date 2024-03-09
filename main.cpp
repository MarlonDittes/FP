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

    g->printGraph();

    long crossing_count = g->countCrossingsMarlon();
    std::cout << "number of crossings in default g: " << crossing_count << std::endl;

    g->setNoPartitioning();

    g->printGraph();
    
    //Testing partitioning
    auto nodes = g->getGraph();
    std::vector<Node*> partition(0);
    partition.push_back(&nodes[0]);
    partition.push_back(&nodes[1]);
    partition.push_back(&nodes[2]);
    partition.push_back(&nodes[5]);
    partition.push_back(&nodes[6]);

    std::vector<Node*> partition2(0);
    partition2.push_back(&nodes[2]);
    partition2.push_back(&nodes[3]);
    partition2.push_back(&nodes[4]);
    partition2.push_back(&nodes[7]);
    partition2.push_back(&nodes[8]);
    
    auto partGraph = createGraphByPartition(g, partition);
    partGraph->printGraph();
    auto nodes1 = partGraph->getGraph();
    for (auto node : nodes1){
        std::cout << node.id << ", previously: " << node.old_id << std::endl;
    }

    auto partGraph2 = createGraphByPartition(g, partition2);
    partGraph2->printGraph();
    auto nodes2 = partGraph2->getGraph();
    for (auto node : nodes2){
        std::cout << node.id << ", previously: " << node.old_id << std::endl;
    }




    /*
    if (g->verifier(verifier)) {
        std::cout << "Graph is valid"<<std::endl;
    }
    else {
        std::cout << "Graph is NOT valid" << std::endl;
    }*/

    // Testing on test set
    /*
    for (int mode = 1; mode<=MODE_COUNT; mode++){
        calculatePerformance("../tiny_test_set", "../performance_data/" + std::to_string(mode), mode);
        std::cout << "Mode " << mode << " done." << std::endl;
    } 
    */

    //readWeightedHyperGraph("../hypergraphs/4000_4_40_5.4.hgr", "../bipartitegraphs/4000_4_40_5.4.gr");

}

