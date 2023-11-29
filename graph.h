#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>

struct Node {
    std::vector<int> X; //name needs to be changed
    int order;
    int id;
    int median = 0;
    std::vector<int> partition;
    int visible_nodes = 0; //offset to visible nodes in X
};

typedef std::vector<std::pair<Node*, Node*>> Twins;

class Graph {
private:
    std::vector<Node> Y;
    std::vector<Node> fixed_nodes;
    std::vector<Node*> Order;
    int size_X;
    int size_Y;
    int edge_number;

public:
    Graph(int size_X, int size_Y, int edge_number);
    void addEdge(int y, int x);
    Node* getNodeByOrder(int order);
    int getSizeOfOrder() { return Order.size(); }
    void printGraph();
    int countCrossings();
    int countCrossingsForPair(int a, int b);
    void sortYArray();
    void swapNodes(int node0, int node1);
    void makeNodeInvisible(int order);
    std::pair<std::vector<Node*>, int> Greedy();
    
    std::pair<int,int> DFS_partition(int start_node_fixed_id, int partition, std::vector<bool>& visited);
    void Partition();
    void APUtil(int start_node_id, std::vector<int>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<int>& isAP);
    void AP();
    void Median_Heuristic();

    bool verifier(Graph check);

    std::vector<Node*> getOrder();
    void setOrder(std::vector<Node*> order);

    
};

Graph* readGraph(std::string graph_file);

std::pair<std::vector<Node*>, int> bruteForce(Graph g);
void outputOrder(std::vector<Node*> order, std::string output);

//reduction:
Twins findTwins();
void cheapReduction();

#endif //GRAPH_H
