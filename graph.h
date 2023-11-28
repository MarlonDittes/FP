#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>

struct Node {
    std::vector<int> X;
    int order;
    int id;
    int median = 0;
    int partition = -1;

    Node() {
        X = std::vector<int>();
    }
};

typedef std::vector<std::pair<Node*, Node*>> Twins;

class Graph {
    private:
    std::vector<Node> Y; 
    std::vector<Node*> Order;
    int size_X;
    int size_Y;
    int edge_number;

    public:
    Graph(int size_X, int size_Y, int edge_number);
    void addEdge(int y, int x);
    Node* getNode(int i);
    void printGraph();
    int countCrossings();
    int countCrossingsForPair(int a, int b);
    void sortYArray();
    void swapNodes(int node0, int node1);
    std::pair<std::vector<Node*>, int> Greedy();
    void DFS_partition();

    void Median_Heuristic();

    bool verifier(Graph check);

    std::vector<Node*> getOrder();
    void setOrder(std::vector<Node*> order);

    //reductions:
    Twins findTwins();
    void cheapReduction();
};

Graph* readGraph(std::string graph_file);

std::pair<std::vector<Node*>, int> bruteForce(Graph g);
void outputOrder (std::vector<Node*> order, std::string output);
#endif //GRAPH_H