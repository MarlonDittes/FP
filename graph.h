#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include "reductions.h"

class general_reduction;

struct Edge {
    int neighbour_id;
    int edge_weight;
};

struct Node {
    // Edge stuff
    std::vector<Edge> edges;
    std::vector<int> neighbors_partition;
    int offset_visible_nodes = 0; //offset to visible nodes in neighbours

    // ID stuff 
    int id;
    int old_id = -1; // refers to previous ID of this node in bigger graph
    int order = -1;
    int median = 0;

    // Partition stuff
    bool isAP = false;
    std::vector<int> partition;


    int multiplier = 1; //REMOVE WHEN EDGE SYSTEM READY
};

struct Partition_Intervall {

    Partition_Intervall() {};

    Partition_Intervall(int low, int high) : interval_high(high), interval_low(low) {};

    std::vector<int> nodes;
    int interval_high = 0;
    int interval_low = 0;
    bool ignore = false;
};

class Graph {
private:
    std::vector<Node> graph;            // Adjacency List
    std::vector<Node*> order_nodes;     // Order of moveable nodes
    int offset_visible_order_nodes = 0; // Offset for accessing visible moveable nodes
    std::vector<std::vector<Node*>> partitions;

    int n0; //size fixed nodes
    int n1; //size movable nodes
    int m; //number of edges
    int activeEdges;

    bool optimal = 0; //for reduction of partitions -> true if partition is optimal

public:

    Graph(int n0, int n1, int m);   //Constructor
    void addEdge(int y, int x);

    // Accessing node information
    Node* getNodeByOrder(int order) { return order_nodes[order]; };
    int getOrderByNode(int node_id) { return graph[node_id].order; }
    void setOrderByNode(int node_id, int newOrder) { graph[node_id].order = newOrder; };
    void setOldID(int node_id, int old_id) { this->graph[node_id].old_id = old_id; };

    // Vector Stuff
    int getSizeOfOrder() { return order_nodes.size(); }
    std::vector<Node*> getOrderNodes() { return this->order_nodes; };
    void setOrderNodes(std::vector<Node*> order);

    std::vector<std::vector<Node*>>& getPartitions() { return this->partitions; };
    std::vector<Node>& getGraph() { return this->graph; };

    // Accessing graph information
    int getN0() { return this->n0; };
    int getN1() { return this->n1; };
    int getM() { return this->m; };
    int getOffsetVisibleOrderNodes() { return this->offset_visible_order_nodes; };
    void setActiveEdges(int count) { this->activeEdges = count; };
    int getActiveEdges() { return this->activeEdges; };

    bool getOptimal() { return this->optimal; };
    void setOptimalTrue() { this->optimal = true; };

    // Printing 
    void printGraph();
    void printGraphByPartitions();

    // Sort
    void sortOrderNodesByOrder();
    void sortNeighbours();
    static bool compareNeighbours(const Edge& a, const Edge& b) { return a.neighbour_id < b.neighbour_id; }

    // Count crossings
    long countCrossingsMarlon();
    long countCrossingsBranching();
    long countCrossings();
    int countCrossingsForPair(int order_a, int order_b);
    long countCrossingsWithEdgeWeights();

    // Swap
    void swapNodes(int order_a, int order_b);
    void swapNodesBranching(int order_a, int order_b);

    // Invisible Node Stuff
    void makeNodeInvisibleMarlon(int order_of_node);
    void makeNodeInvisible(int order);
    void makeNodeInvisibleBranching(int order);

    void makeNodeVisibleMarlon();
    void makeNodeVisible(int order_of_node);

    // Partitioning
    void setNoPartitioning();

    void Interval_Partitioning(int start_node_id, int end_node_id);
    std::pair<int, int> DFSforPartition(int start_node_fixed_id, int partition, std::vector<bool>& visited);
    void Partition();
    void APUtil(int start_node_id, std::vector<bool>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<bool>& isAP);
    void APUtilIterative(int start_node_id, std::vector<bool>& visited, std::vector<int>& disc, std::vector<int>& low, std::vector<int>& parent, std::vector<bool>& isAP);
    int CreatePartitionsVector(int start_node_fixed_id, int& partition_id, std::vector<bool>& visited);
    void AP();
    void AP_Intervall();
    int CreatePartitionsVectorNew(int start_node_fixed_id, int& partition_id, std::vector<bool>& visited, std::vector<std::vector<Node*>>& partitions_temp);
    void DFS_AP_nodes(int start_node, std::vector<bool>& visited, std::vector<Partition_Intervall>& partition_intervall);

    // Solving
    std::pair<std::vector<Node*>, long> Greedy();
    void MedianHeuristicMarlon();
    void Median_Heuristic();
    void Sorted_straight_line_reduction();
    bool DFS_for_sorted_straight_line(int start_node, std::vector<bool>& visited);

    // Verifier
    bool verifier(Graph check);
};

// Brute Force
std::pair<std::vector<Node*>, long> bruteForce(Graph* g);
std::pair<std::vector<Node*>, long> bruteForceOnSubgraph(Graph* g, int begin, int end);

// Branching Stuff
Graph* createGraphByPartition(Graph* g, std::vector<Node*> partition);
std::pair<std::vector<Node*>, long> branching(Graph* g, std::vector<general_reduction*> reductionTypes);
std::pair<std::vector<Node*>, long> BranchAndReduce(Graph* g, std::vector<general_reduction*> reductionTypes);

// OLD SHAI BRANCHING -> REMOVE?
void Branch_and_Bound(Graph* G);
void exploreBranch(Graph G_branch, Graph& G_original, int depth, int& best_solution, std::vector<Node*>& best_configuration);

#endif //GRAPH_H

