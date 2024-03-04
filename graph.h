#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include "stack.h"

struct Node {
    std::vector<int> neighbours;
    int order = -1;
    Stack ID_Stack;
    int median = 0;
    int partition = -1;
    //std::vector<int> partition;
    int offset_visible_nodes = 0; //offset to visible nodes in neighbours
    //std::vector<int> offset_neighbours_partition;
};

typedef std::vector<std::pair<Node*, Node*>> Twins;

class Graph {
private:
    std::vector<Node*> graph;
    std::vector<Node*> order_nodes;
    std::vector<std::vector<Node*>> partitions;
    int n0; //size fixed nodes
    int n1; //size movable nodes
    int m; //number of edges
    //std::vector<std::vector<int>> partitions;
    std::vector<Node*> activeNodes;

public:
    Graph(int n0, int n1, int m);
    void addEdge(int y, int x);

    Node* getNodeByOrder(int order) { return order_nodes[order]; };
    int getOrderByNode(int node_id) { return graph[node_id]->order; }
    int getSizeOfOrder() { return order_nodes.size(); }
    std::vector<Node*> getOrderNodes() { return this->order_nodes; };
    void setOrderNodes(std::vector<Node*> order) { this->order_nodes = order; };
    // New getters, do we need?
    std::vector<std::vector<Node*>> getPartitions() { return this->partitions; };
    int getN0() { return this->n0; };
    int getN1() { return this->n1; };
    int getM() { return this->m; };
    std::vector<Node*> getGraph() { return this->graph; };
    std::vector<Node*> getActiveNodes() { return this->activeNodes;}

    void printGraph();
    long countCrossings();
    long countCrossingsBranching();
    int countCrossingsForPair(int order_a, int order_b);
    void sortNeighbours();
    void swapNodes(int order_a, int order_b);
    void swapNodesBranching(int order_a, int order_b);

    void makeNodeInvisible(int order);
    void makeNodeInvisibleBranching(int order);
    void Interval_Partitioning(int start_node_id, int end_node_id);

    std::pair<std::vector<Node*>, long> Greedy();
    std::pair<int, int> DFSforPartition(int start_node_fixed_id, int partition, std::vector<bool>& visited);
    void Partition();
    void APUtil(int start_node_id, std::vector<int>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<int>& isAP);
    void AP();
    void Median_Heuristic();
    void Sorted_straight_line_reduction();
    void makeNodeVisible(int order_of_node);
    bool DFS_for_sorted_straight_line(int start_node, std::vector<bool>& visited);
    bool verifier(Graph check);
};

std::pair<std::vector<Node*>, long> bruteForce(Graph* g);
std::pair<std::vector<Node*>, long> bruteForceOnSubgraph(Graph* g, int begin, int end);
// New stuff
std::pair<std::vector<Node*>, long> BranchAndReduce(Graph G);
bool reduceCliqueReduction(Graph* g);

void Branch_and_Bound(Graph* G);
void exploreBranch(Graph G_branch, Graph& G_original, int depth, int& best_solution, std::vector<Node*>& best_configuration);

//reduction:
Twins findTwins();
void cheapReduction();

#endif //GRAPH_H
