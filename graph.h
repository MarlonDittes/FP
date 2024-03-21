#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include "reductions.h"

class general_reduction;

struct Node {
    std::vector<int> neighbours;
    std::vector<int> neighbors_partition;
    int order = -1;
    int id;
    int median = 0;
    bool isAP = false;
    //int partition = -1;
    std::vector<int> partition;
    int offset_visible_nodes = 0; //offset to visible nodes in neighbours
    //std::vector<int> offset_neighbours_partition;
    int old_id = -1;
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
    //std::vector<std::vector<int>> partitions;
    int activeEdges;
    int offset_visible_order_nodes = 0;
    bool optimal = 0; //for reduction of partitions -> true if partition is optimal

public:
    int global_partition = 0;
    Graph(int n0, int n1, int m);
    void addEdge(int y, int x);

    Node* getNodeByOrder(int order) { return order_nodes[order]; };
    int getOrderByNode(int node_id) { return graph[node_id].order; }
    void setOrderByNode(int node_id, int newOrder) { graph[node_id].order = newOrder; };
    int getSizeOfOrder() { return order_nodes.size(); }
    std::vector<Node*> getOrderNodes() { return this->order_nodes; };
    void setOrderNodes(std::vector<Node*> order);
    // New getters, do we need?
    std::vector<std::vector<Node*>>& getPartitions() { return this->partitions; };
    int getN0() { return this->n0; };
    int getN1() { return this->n1; };
    int getM() { return this->m; };
    int getOffsetVisibleOrderNodes() { return this->offset_visible_order_nodes; };
    std::vector<Node>& getGraph() { return this->graph; };
    void setOldID(int node_id, int old_id) { this->graph[node_id].old_id = old_id; };
    void setActiveEdges (int count) { this->activeEdges = count;};
    int getActiveEdges () { return this->activeEdges;};

    void sortOrderNodesByOrder();

    long countCrossingsMarlon();
    bool getOptimal() { return this->optimal; };
    void setOptimalTrue() { this->optimal = true; };

    void printGraph();
    void printGraphByPartitions();

    long countCrossings();
    long countCrossingsBranching();
    int countCrossingsForPair(int order_a, int order_b);
    long countCrossingsWithMultiplier();
    void sortNeighbours();
    void swapNodes(int order_a, int order_b);
    void swapNodesBranching(int order_a, int order_b);

    void makeNodeInvisibleMarlon(int order_of_node);
    void makeNodeInvisible(int order);
    void makeNodeInvisibleBranching(int order);
    void makeNodeVisibleMarlon();
    void makeNodeVisible(int order_of_node);
    void Interval_Partitioning(int start_node_id, int end_node_id);

    std::pair<std::vector<Node*>, long> Greedy();
    std::pair<int, int> DFSforPartition(int start_node_fixed_id, int partition, std::vector<bool>& visited);
    void Partition();
    void APUtil(int start_node_id, std::vector<bool>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<bool>& isAP);
    void APUtilIterative(int start_node_id, std::vector<bool>& visited, std::vector<int>& disc, std::vector<int>& low, std::vector<int>& parent, std::vector<bool>& isAP);
    int CreatePartitionsVector(int start_node_fixed_id, int& partition_id, std::vector<bool>& visited);
    void AP();
    void MedianHeuristicMarlon();
    void Median_Heuristic();
    void Sorted_straight_line_reduction();
    bool DFS_for_sorted_straight_line(int start_node, std::vector<bool>& visited);
    bool verifier(Graph check);

    void setNoPartitioning();
};

std::pair<std::vector<Node*>, long> bruteForce(Graph* g);
std::pair<std::vector<Node*>, long> bruteForceOnSubgraph(Graph* g, int begin, int end);
// New stuff
Graph* createGraphByPartition(Graph* g, std::vector<Node*> partition);
std::pair<std::vector<Node*>, long> branching(Graph* g, std::vector<general_reduction*> reductionTypes);
std::pair<std::vector<Node*>, long> BranchAndReduce(Graph* g, std::vector<general_reduction*> reductionTypes);

void Branch_and_Bound(Graph* G);
void exploreBranch(Graph G_branch, Graph& G_original, int depth, int& best_solution, std::vector<Node*>& best_configuration);



#endif //GRAPH_H

