#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>

struct Node {
    std::vector<int> X;
    int order;
    int id;

    Node() {
        X = std::vector<int>();
    }
};

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
    void sortYArray();

    std::vector<Node*> getOrder();
    void setOrder(std::vector<Node*> order);
};

Graph* readGraph(std::string graph_file);

std::pair<std::vector<Node*>, int> bruteForce(Graph g);
void outputOrder (std::vector<Node*> order, std::string output);
#endif //GRAPH_H