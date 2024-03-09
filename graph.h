#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>

struct Node {
    std::vector<int> neighbours;
    int order = -1;
    int id = -1;
    int median = 0;
    int partition = -1;
    int offset_visible_nodes = 0; //offset to visible nodes in neighbours
    int multiplier = 1;
};

class Graph {
private:
    std::vector<Node> graph;
    std::vector<Node*> order_nodes;
    std::vector<std::vector<Node*>> partitions;
    int n0; //size fixed nodes
    int n1; //size movable nodes
    int m; //number of edges
    bool optimal = 0; //for reduction of partitions -> true if partition is optimal

public:
    Graph(int n0, int n1, int m);
    void addEdge(int y, int x);

    Node* getNodeByOrder(int order) { return order_nodes[order]; };
    int getSizeOfOrder() { return order_nodes.size(); }
    std::vector<Node*> getOrderNodes() { return this->order_nodes; };
    void setOrderNodes(std::vector<Node*> order) { this->order_nodes = order; };
    // New getters, do we need?
    std::vector<std::vector<Node*>> getPartitions() { return this->partitions; };
    int getN0() { return this->n0; };
    int getN1() { return this->n1; };
    int getM() { return this->m; };
    std::vector<Node> getGraph() { return this->graph; };
    bool getOptimal() { return this->optimal; };
    void setOptimalTrue() { this->optimal = true; };

    void printGraph();
    void printGraphByPartitions();
    long countCrossings();
    int countCrossingsForPair(int order_a, int order_b);
    long countCrossingsWithMultiplier();
    void sortNeighbours();
    void swapNodes(int order_a, int order_b);
    void makeNodeInvisible(int order);

    std::pair<std::vector<Node*>, long> Greedy();
    void Median_Heuristic();
    bool verifier(Graph check);
    std::pair<int,int> DFSforPartition(int start_node_fixed_id, int partition, std::vector<bool>& visited);
    void Partition();
    void APUtil(int start_node_id, std::vector<int>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<int>& isAP);
    void AP();
    bool DFS_for_sorted_straight_line(int start_node, std::vector<bool>& visited);
    void Sorted_straight_line_reduction();
};

std::pair<std::vector<Node*>, long> bruteForce(Graph* g);
std::pair<std::vector<Node*>, long> bruteForceOnSubgraph(Graph* g, int begin, int end);
Graph* createGraphByPartition(Graph* g, std::vector<Node*> partition);
// New stuff
std::pair<std::vector<Node*>, long> BranchAndReduce(Graph G);
bool reduceCliqueReduction(Graph* g);

#endif //GRAPH_H
