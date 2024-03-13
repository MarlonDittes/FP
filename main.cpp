#include <iostream>
#include "graph.h"
#include "io.h"
#include "stack.h"
#include "performance.h"

int main(int argc, char* argv[]) {
    std::string graph_file = argv[1];
    std::cout << graph_file << std::endl;

    
    Graph* g = readGraph(graph_file);
    Graph verifier = *g;

    g->printGraph();

    long crossing_count = g->countCrossingsMarlon();
    std::cout << "number of crossings in default g: " << crossing_count << std::endl;

    std::vector<Reduction> reductions;
    

    /*
    //Testing partitioning

    //g->makeNodeInvisibleMarlon(0);
    
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
    auto order = branching(partGraph, reductions);

    auto partGraph2 = createGraphByPartition(g, partition2);
    auto order2 = branching(partGraph2, reductions);

    std::vector<Node*> solution(0);
    solution.insert(solution.end(), order.first.begin(), order.first.end());
    solution.insert(solution.end(), order2.first.begin(), order2.first.end());


    std::vector<Node*> oldOrder = g->getOrderNodes();
    nodes = g->getGraph();
    std::vector<Node*> newOrder(oldOrder.begin(), oldOrder.begin() + g->getOffsetVisibleOrderNodes());
    for (auto node : solution){
        newOrder.push_back(&nodes[node->old_id]);
    }
    g->setOrderNodes(newOrder);

    g->printGraph();
    */

    
    //Testing Branching
    
    auto result = BranchAndReduce(g, reductions);
    g->setOrderNodes(result.first);

    g->printGraph();
    




    /*
    if (g->verifier(verifier)) {
        std::cout << "Graph is valid"<<std::endl;
    }
    else {
        std::cout << "Graph is NOT valid" << std::endl;
    }*/

    //Testing on test set

    std::vector<int> modes;
    modes.push_back(1);
    modes.push_back(2);
    modes.push_back(3);
    modes.push_back(4);
    modes.push_back(5);

    /*
    for (auto mode : modes){
        calculatePerformance("../tiny_test_set", "../performance_data/" + std::to_string(mode), mode);
        std::cout << "Mode " << mode << " done." << std::endl;
    }
    */ 
}

