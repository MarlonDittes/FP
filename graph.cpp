#include "graph.h"
#include "stack.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <utility>
#include <cmath>


//Constructor dependent on size of Y
Graph::Graph(int size_X, int size_Y, int edge_number) : size_X(size_X), size_Y(size_Y), edge_number(edge_number), fixed_nodes(size_X){
    std::cout << "size_X: " << size_X << std::endl;
    std::cout << "assert" << std::endl;
    Y = std::vector<Node>(size_Y);
    Order = std::vector<Node*>(size_Y);

    //initialize id and order in the Y and Order Array
    for (int i = 0; i < Order.size(); i++) {
        Y[i].id = size_X + i + 1; //set id to start at n0+1 to n0+n1;
        Y[i].order = i;
        //Order[i] = Y[i].id //link order to id 
        Order[i] = &Y[i]; //link order to index (not id!)
    }

    for (int i = 0; i < this->size_X; i++) {
        fixed_nodes[i].id = i + 1;
    }
}

//adding Edges from Y to X
void Graph::addEdge(int x, int y) {
    //!offset of X because y id starts at n0+1 -> this function works only for input data
    Y[y - size_X - 1].X.push_back(x); // in vector Y add adjacency (y,x) with x of X
    fixed_nodes[x - 1].X.push_back(y);
}

Node* Graph::getNodeByOrder(int order) {
    return Order[order];
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

    for (int i = 0; i < Y.size(); i++) {
        std::cout << "Node Id: " << Y[i].id << " Partition : ";
        for (int j = 0; j < Y[i].partition.size(); j++) {
            std::cout << Y[i].partition[j] << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < fixed_nodes.size(); i++) {
        std::cout << "Node Id: " << fixed_nodes[i].id << " Partition : ";
        for (int j = 0; j < fixed_nodes[i].partition.size(); j++) {
            std::cout<<fixed_nodes[i].partition[j] << " ";
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
    for (int u = 0; u < Order.size() - 1; u++) {
        //iterate through Y+1
        for (int v = u + 1; v < Order.size(); v++) {
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

int Graph::countCrossingsForPair(int a, int b) { //for a < b
    int crossings = 0;
    for (auto& i : Order[a]->X) {
        for (auto& j : Order[b]->X) {
            if (i > j) {
                crossings++;
            }
        }
    }
    return crossings;
}

int Graph::countCrossingsForPair(int a, int b) { //for a < b
    int crossings = 0;
    for (auto& i : Order[a]->X) {
        for (auto& j : Order[b]->X) {
            if (i > j) {
                crossings++;
            }
        }
    }
    return crossings;
}

std::vector<Node*> Graph::getOrder() {
    return this->Order;
}

void Graph::setOrder(std::vector<Node*> order) {
    this->Order = order;
}

void Graph::swapNodes(int node0, int node1) {
    std::swap(Order[node0]->order, Order[node1]->order);
    std::swap(Order[node0], Order[node1]);
}

void Graph::makeNodeInvisible(int order_of_node) {
    Node* node = Order[order_of_node];
    for (int i = node->visible_nodes; i < node->X.size(); i++) {
        Node* neighbour = &fixed_nodes[node->X[i]];

        //make node invisible in all adjacencie lists of fixed neighbours
        for (int j = neighbour->visible_nodes; j < neighbour->X.size(); j++) {
            if (neighbour->X[j] != node->id) {
                continue;
            }
            //add node to invisible nodes and add visible counter
            std::swap(neighbour->X[j], neighbour->X[neighbour->visible_nodes]);
            neighbour->visible_nodes++;
        }
    }

    //set visible counter to end of adjacency array -> as if all edges removed
    node->visible_nodes = node->X.size();
}

std::pair<std::vector<Node*>, int> Graph::Greedy() {
    int crossings;
    for (int i = 0; i < Order.size() - 1; i++) {
        int crossings_old = this->countCrossings();
        this->swapNodes(i, i + 1);
        int crossings_new = this->countCrossings();
        if (crossings_old < crossings_new) {
            this->swapNodes(i, i + 1);
        }
        if (i == Order.size() - 2) {
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

std::pair<int, int> Graph::DFS_partition(int start_node_fixed_id, int partition ,std::vector<bool>& visited) {    
    // need to make sure at the beginning that all nodes with degree 0 are removed. we can do this since by setting their
    // order to be the right most position and adding an ignore marker.
    Stack stack;

    stack.push(start_node_fixed_id);

    int start_node_interval = start_node_fixed_id;
    int last_node_interval = -1;

    //move this to a seperate function called DFS
    while (!stack.isEmpty()) {
        int current_node = stack.peek();
        stack.pop();

        if (current_node < size_X && current_node > last_node_interval) {
            last_node_interval = current_node;
        }

        if (!visited[current_node]) {
            std::cout << "Visited " << current_node << " ";
            visited[current_node] = true;
            if (current_node >= size_X) {
                // Y[current_node - size_X].partition = partition;
                Y[current_node - size_X].partition.push_back(partition);
            }
            else {
                fixed_nodes[current_node].partition.push_back(partition);
            }
        }

        if (current_node >= size_X) {
            for (int i = 0; i < Y[current_node - size_X].X.size(); i++) {
                if (!visited[Y[current_node - size_X].X[i] - 1]) {
                    stack.push(Y[current_node - size_X].X[i] - 1);
                }
            }
        }
        else {
            for (int i = 0; i < fixed_nodes[current_node].X.size(); i++) {
                if (!visited[fixed_nodes[current_node].X[i] - 1]) {
                    stack.push(fixed_nodes[current_node].X[i] - 1);
                }
            }
        }
    }

    return std::make_pair(start_node_interval, last_node_interval);
}


void Graph::Partition()
{
    int start_node_fix = fixed_nodes[0].id - 1;
    int partition = 0;
    int rounds_counter = 0;

    std::vector<bool> visited(this->size_X + this->size_Y, false);

    while (rounds_counter < this->size_X) {
        
        std::pair<int, int> end_start_node_interval = DFS_partition(start_node_fix, partition, visited);
        rounds_counter++;
        int next_fixed_node = end_start_node_interval.first + 1;
        
        while (next_fixed_node <= end_start_node_interval.second) {
            if (visited[next_fixed_node]) {
                next_fixed_node++;
                rounds_counter++;
            }
            else {
                std::pair<int,int> end_start_node_interval_2 = DFS_partition(next_fixed_node, partition, visited);
                rounds_counter++;
                next_fixed_node = end_start_node_interval_2.first + 1;
            }
        }
        start_node_fix = next_fixed_node;
        partition++;
    }
}


// A recursive function that find articulation 
// points using DFS traversal
// adj[] --> Adjacency List representation of the graph
// start_node_id --> The vertex to be visited next
// visited[] --> keeps track of visited vertices
// disc[] --> Stores discovery times of visited vertices
// low[] -- >> earliest visited vertex (the vertex with minimum
// discovery time) that can be reached from subtree
// rooted with current vertex
// parent --> Stores the parent vertex in DFS tree
// isAP[] --> Stores articulation points
void Graph::APUtil(int start_node_id, std::vector<int>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<int>& isAP)
{
    // Count of children in DFS Tree
    int children = 0;

    // Mark the current node as visited
    visited[start_node_id] = true;

    // Initialize discovery time and low value
    disc[start_node_id] = low[start_node_id] = ++time;

    if (start_node_id >= size_X) {
        // Go through all vertices adjacent to this
        for (int i = 0; i < Y[start_node_id - size_X].X.size(); i++) {
            // If v is not visited yet, then make it a child of u
            // in DFS tree and recur for it
            if (!visited[Y[start_node_id - size_X].X[i] - 1]) {
                children++;
                APUtil(Y[start_node_id - size_X].X[i] - 1, visited, disc, low, time, start_node_id, isAP);

                // Check if the subtree rooted with v has
                // a connection to one of the ancestors of u
                low[start_node_id] = std::min(low[start_node_id], low[Y[start_node_id - size_X].X[i] - 1]);

                // If start_node_id is not root and low value of one of
                // its child is more than discovery value of u.
                if (parent != -1 && low[Y[start_node_id - size_X].X[i] - 1] >= disc[start_node_id])
                    isAP[start_node_id] = true;
            }

            // Update low value of start_node_id for parent function calls.
            else if (Y[start_node_id - size_X].X[i] - 1 != parent)
                low[start_node_id] = std::min(low[start_node_id], disc[Y[start_node_id - size_X].X[i] - 1]);
        }
    }
    else {
        // Go through all vertices adjacent to this
        for (int i = 0; i < fixed_nodes[start_node_id].X.size(); i++) {
            // If v is not visited yet, then make it a child of u
            // in DFS tree and recur for it
            if (!visited[fixed_nodes[start_node_id].X[i] - 1]) {
                children++;
                APUtil(fixed_nodes[start_node_id].X[i] - 1, visited, disc, low, time, start_node_id, isAP);

                // Check if the subtree rooted with v has
                // a connection to one of the ancestors of u
                low[start_node_id] = std::min(low[start_node_id], low[fixed_nodes[start_node_id].X[i] - 1]);

                // If start_node_id is not root and low value of one of
                // its child is more than discovery value of u.
                if (parent != -1 && low[fixed_nodes[start_node_id].X[i] - 1] >= disc[start_node_id])
                    isAP[start_node_id] = true;
            }

            // Update low value of start_node_id for parent function calls.
            else if (fixed_nodes[start_node_id].X[i] - 1 != parent)
                low[start_node_id] = std::min(low[start_node_id], disc[fixed_nodes[start_node_id].X[i] - 1]);
        }
    }

    // If start_node_id is root of DFS tree and has two or more children.
    if (parent == -1 && children > 1)
        isAP[start_node_id] = true;
}

void Graph::AP()
{
    std::vector<int> disc(this->size_X + this->size_Y, 0);
    std::vector<int> low(this->size_X + this->size_Y);
    std::vector<int> visited(this->size_X + this->size_Y, false);
    std::vector<int> isAP(this->size_X + this->size_Y, false);
    int time = 0, par = -1;

    // Adding this loop so that the
    // code works even if we are given
    // disconnected graph
    for (int start_node_id= 0; start_node_id < this->size_X + this->size_Y; start_node_id++)
        if (!visited[start_node_id]) {
            APUtil(start_node_id, visited, disc, low, time, par, isAP);
        }

    // Printing the APs
    for (int start_node_id = 0; start_node_id < this->size_X + this->size_Y; start_node_id++) {
        if (isAP[start_node_id]) {
            std::cout << "Node id: " << start_node_id + 1 << " is AP" << std::endl;
        }
    }

    int partition_offset = 0;
    for (int i = 0; i < this->size_X; i++) {
        for (int j = 0; j < fixed_nodes[i].partition[j]; j++) {
            fixed_nodes[i].partition[j] = fixed_nodes[i].partition[j] + partition_offset;

        }
        if (isAP[i]) {
            partition_offset++;


        }
    }

}

//typedef std::vector<std::pair<Node*, Node*>> Twins;

Twins findTwins(Graph* g) {
    Twins twins;
    for (int i = 0; i < g->getSizeOfOrder()-1; i++) {
        for (int j = i+1; j < g->getSizeOfOrder(); j++) {
            if (g->getNodeByOrder(i)->X.size() != g->getNodeByOrder(j)->X.size()) {
                continue;
            }

            for (int k = 0; k < g->getNodeByOrder(i)->X.size(); k++) {
                if (g->getNodeByOrder(i)->X[k] != g->getNodeByOrder(j)->X[k]) {
                    break;
                }

                //found twins
                if (k == g->getNodeByOrder(i)->X.size()-1) {
                    std::pair<Node*, Node*> twin = std::make_pair(g->getNodeByOrder(i), g->getNodeByOrder(j)); 
                    twins.push_back(twin); //twin_reduce()
                }
            }
        }
    }
    return twins;
}

/*
void Graph::cheapReduction() {
    for (int a = 0; a < Order.size()-1; a++) {
        for (int b = a+1; b < Order.size(); b++) {
            int a_smaller_b = countCrossingsForPair(a, b);
            int b_smaller_a = countCrossingsForPair(b, a);
            if (a_smaller_b == 0 || b_smaller_a == 0) {

            }
        }
    }
}
*/

bool compareNodePointers(Node* a, Node* b) {
    return a->order < b->order;
}

std::pair<std::vector<Node*>, int> bruteForce(Graph g) {
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

        if (crossings < bestCrossings) {
            bestOrder = baseOrder;
            bestCrossings = crossings;
        }

    } while (std::next_permutation(baseOrder.begin(), baseOrder.end(), compareNodePointers));

    return std::make_pair(bestOrder, bestCrossings);
}