#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>

struct Node {
    std::vector<int> neighbours; //name needs to be changed
    int order = -1;
    int id = -1;
    int median = 0;
    std::vector<int> partition;
    int offset_visible_nodes = 0; //offset to visible nodes in neighbours
};

typedef std::vector<std::pair<Node*, Node*>> Twins;

class Graph {
private:
    std::vector<Node> graph;
    std::vector<Node*> order_nodes;
    int n0; //size fixed nodes
    int n1; //size movable nodes
    int m; //number of edges

public:
    Graph(int n0, int n1, int m);
    void addEdge(int y, int x);

    Node* getNodeByOrder(int order) { return order_nodes[order]; };
    int getSizeOfOrder() { return order_nodes.size(); }
    std::vector<Node*> getOrderNodes() { return this->order_nodes; };
    void setOrderNodes(std::vector<Node*> order) { this->order_nodes = order; } ;

    void printGraph();
    long countCrossings();
    int countCrossingsForPair(int order_a, int order_b);
    void sortNeighbours();
    void swapNodes(int order_a, int order_b);
    void makeNodeInvisible(int order);

    std::pair<std::vector<Node*>, long> Greedy();
    std::pair<int,int> DFSforPartition(int start_node_fixed_id, int partition, std::vector<bool>& visited);
    void Partition();
    void APUtil(int start_node_id, std::vector<int>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<int>& isAP);
    void AP();
    void Median_Heuristic();
    void Sorted_straight_line_reduction();
    bool DFS_for_sorted_straight_line(int start_node, std::vector<bool>& visited);
    bool verifier(Graph check);
};

std::pair<std::vector<Node*>, long> bruteForce(Graph g);
long factorial(int n);
std::pair<std::vector<Node*>, long> bruteForceParallel(Graph g);

//reduction:
Twins findTwins();
void cheapReduction();

#endif //GRAPH_H
