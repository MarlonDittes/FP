#include "graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <utility>
#include <cmath>

//Constructor dependent on size of Y
Graph::Graph(int size_X, int size_Y, int edge_number) : size_X(size_X), size_Y(size_Y), edge_number(edge_number) {
    std::cout << "size_X: " << size_X << std::endl;
    std::cout << "assert" << std::endl;
    Y = std::vector<Node>(size_Y);
    Order = std::vector<Node*>(size_Y);

    //initialize id and order in the Y and Order Array
    for (int i = 0; i < Order.size(); i++) {
        Y[i].id = size_X+i+1; //set id to start at n0+1 to n0+n1;
        Y[i].order = i; 
        //Order[i] = Y[i].id //link order to id 
        Order[i] = &Y[i]; //link order to index (not id!)
    }
}

//adding Edges from Y to X
void Graph::addEdge(int x, int y) {
    //!offset of X because y id starts at n0+1 -> this function works only for input data
    Y[y-size_X-1].X.push_back(x); // in vector Y add adjacency (y,x) with x of X 
}

Node* Graph::getNode(int i) {
    return &Y[i];
}

void Graph::printGraph() {

    //print in Order of Y (B)
    for (int i = 0; i < Order.size(); i++) {
        std::cout << Order[i]->id << " -> ";
        for (int j = 0; j < Order[i]->X.size(); j++) {
            std::cout << Order[i]->X[j] << " ";
        }
        std::cout << std::endl;
    }
}

void Graph::sortYArray() {
    for (int i = 0; i < Y.size(); i++) {
        std::sort(Y[i].X.begin(), Y[i].X.end());
    }
}

int Graph::countCrossings() {
    if (edge_number == 0) {
        return 0;
    }
    int crossings = 0;
    //iterate through Y
    for (int u = 0; u < Order.size()-1; u++) {
        //iterate through Y+1
        for (int v = u+1; v < Order.size(); v++) {
            for (auto& i : Order[u]->X) {
                for (auto& j : Order[v]->X) {
                    if (i > j) {
                        crossings++;
                    }
                }
            }                
        }
    }
    return crossings;
}

std::vector<Node*> Graph::getOrder(){
    return this->Order;
}

void Graph::setOrder(std::vector<Node*> order){
    this->Order = order;
}

void Graph::swapNodes(int node0, int node1) {
    std::swap(Order[node0]->order, Order[node1]->order);
    std::swap(Order[node0], Order[node1]);
}

std::pair<std::vector<Node*>, int> Graph::Greedy() {
    int crossings;          
    for (int i = 0; i < Order.size()-1; i++) {
        int crossings_old = this->countCrossings();
        this->swapNodes(i, i+1);
        int crossings_new = this->countCrossings();
        if (crossings_old < crossings_new) {
            this->swapNodes(i, i+1);
        }
        if (i == Order.size()-2) {
            crossings = std::min(crossings_new, crossings_old);
        }
    }

    return std::make_pair(this->Order, crossings);
}

void Graph::Median_Heuristic()
{
    for (int i = 0; i < Order.size(); i++) {
        for (int j = 0; j < Order[i]->X.size(); j++) {
            Order[i]->median += Order[i]->X[j];
        }

        Order[i]->median = Order[i]->median / Order[i]->X.size();
        std::cout << "id: " << Order[i]->id << " median: " << Order[i]->median << std::endl;
    }
    
    std::sort(Order.begin(), Order.end(), [](const Node* a, const Node* b) {
        return a->median < b->median;
    });
    
    for (int i = 0; i < Order.size(); i++) {
        Order[i]->order = i;
    }
}

bool Graph::verifier(Graph check)
{
    std::vector<bool> unique_id_vec(this->Order.size(), false);

    for (int i = 0; i < Order.size(); i++) {
        
        //check every node id comes exactly once
        if (unique_id_vec[this->Order[i]->id - this->size_X - 1]) { return false; }
        unique_id_vec[this->Order[i]->id - this->size_X - 1] = true;
        
        //check for "illegal node id's"
        if (this->Order[i]->id > (check.size_X + check.size_Y)) { return false; }
        
        //neighbours of each node are the same as the graph at the beginning
        for (int j = 0; j < Order[i]->X.size(); j++) {
            int id_check = this->Order[i]->id;
            if (this->Order[i]->X[j] != check.Y[id_check - this->size_X - 1].X[j]) { return false; }
        }
    }

    //a node id disappeared from the Order array.
    for (int i = 0; i < unique_id_vec.size(); i++) {
        if (!unique_id_vec[i]) { return false; }
    }
    return true;
}

/*void Graph::DFS_partition() {
    int partition = 0;

    std::vector
}*/


bool compareNodePointers(Node* a, Node* b) {
    return a->order < b->order;
}

std::pair<std::vector<Node*>, int> bruteForce(Graph g){
    std::vector<Node*> baseOrder = g.getOrder();
    std::sort(baseOrder.begin(), baseOrder.end(), compareNodePointers);

    std::vector<Node*> bestOrder = baseOrder;
    int bestCrossings = g.countCrossings();

    do {
        Graph tmp = g;
        tmp.setOrder(baseOrder);
        // Display the current permutation
        /*std::cout << "Permutation:";
        for (const Node* node : baseOrder) {
            std::cout << " Order: " << node->order << ", ID: " << node->id;
        }
        std::cout << "\n";
*/
        int crossings = tmp.countCrossings();
        //std::cout << "Crossings: " << crossings << std::endl;

        if (crossings < bestCrossings){
            bestOrder = baseOrder;
            bestCrossings = crossings;
        }

    } while (std::next_permutation(baseOrder.begin(), baseOrder.end(), compareNodePointers));

    return std::make_pair(bestOrder, bestCrossings);
}