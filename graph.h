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
};

typedef std::vector<std::pair<Node*, Node*>> Twins;

class Graph {
private:
    std::vector<Node> graph;
    std::vector<Node*> order_nodes;
    std::vector<std::vector<Node*>> partitions;
    int n0; //size fixed nodes
    int n1; //size movable nodes
    int m; //number of edges
    //std::vector<std::vector<int>> partitions;

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

    void printGraph();
    long countCrossings();
    int countCrossingsForPair(int order_a, int order_b);
    void sortNeighbours();
    void swapNodes(int order_a, int order_b);
    void makeNodeInvisible(int order);

    std::pair<std::vector<Node*>, int> Greedy();
    std::pair<int,int> DFSforPartition(int start_node_fixed_id, int partition, std::vector<bool>& visited);
    void Partition();
    void APUtil(int start_node_id, std::vector<int>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<int>& isAP);
    void AP();
    bool DFS_for_sorted_straight_line(int start_node, std::vector<bool>& visited);
    void Sorted_straight_line_reduction();
};

std::pair<std::vector<Node*>, long> bruteForce(Graph* g);
std::pair<std::vector<Node*>, long> bruteForceOnSubgraph(Graph* g, int begin, int end);
// New stuff
std::pair<std::vector<Node*>, long> BranchAndReduce(Graph G);
bool reduceCliqueReduction(Graph* g);

//reduction:
Twins findTwins();
void cheapReduction();

#endif //GRAPH_H
