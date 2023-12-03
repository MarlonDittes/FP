#include "graph.h"
#include "stack.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <utility>
#include <cmath>
#include <omp.h>


//Constructor dependent on size of movable_nodes
Graph::Graph(int n0, int n1, int m) : n0(n0), n1(n1), m(m), fixed_nodes(n0), movable_nodes(n1), order_nodes(n1) {

    //initialize id and order in the movable_nodes and order_nodes Array
    for (int i = 0; i < order_nodes.size(); i++) {
        movable_nodes[i].id = n0 + i + 1; //set id to start at n0+1 to n0+n1;
        movable_nodes[i].order = i;
        //order_nodes[i] = movable_nodes[i].id //link order to id 
        order_nodes[i] = &movable_nodes[i]; //link order to index (not id!)
    }

    for (int i = 0; i < this->n0; i++) {
        fixed_nodes[i].id = i + 1;
    }
}

//adding Edges from movable_nodes to neighbours
void Graph::addEdge(int x, int y) {
    //!offset of neighbours because y id starts at n0+1 -> this function works only for input data
    movable_nodes[y - n0 - 1].neighbours.push_back(x); // in vector movable_nodes add adjacency (y,x) with x of X
    fixed_nodes[x - 1].neighbours.push_back(y);
}

void Graph::printGraph() {

    //print in order_nodes of movable_nodes (B)
    for (int i = 0; i < order_nodes.size(); i++) {
        std::cout << order_nodes[i]->id << " -> ";
        for (int j = 0; j < order_nodes[i]->neighbours.size(); j++) {
            std::cout << order_nodes[i]->neighbours[j] << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < movable_nodes.size(); i++) {
        std::cout << "Node Id: " << movable_nodes[i].id << " Partition : ";
        for (int j = 0; j < movable_nodes[i].partition.size(); j++) {
            std::cout << movable_nodes[i].partition[j] << " ";
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

void Graph::sortMovableNodes() {
    for (int i = 0; i < movable_nodes.size(); i++) {
        std::sort(movable_nodes[i].neighbours.begin(), movable_nodes[i].neighbours.end());
    }
}

int Graph::countCrossings() {
    if (m == 0) {
        return 0;
    }
    int crossings = 0;
    //iterate through movable_nodes
    for (int u = 0; u < order_nodes.size() - 1; u++) {
        //iterate through movable_nodes+1
        for (int v = u + 1; v < order_nodes.size(); v++) {
            for (auto& i : order_nodes[u]->neighbours) {
                for (auto& j : order_nodes[v]->neighbours) {
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
    for (auto& i : order_nodes[a]->neighbours) {
        for (auto& j : order_nodes[b]->neighbours) {
            if (i > j) {
                crossings++;
            }
        }
    }
    return crossings;
}

void Graph::swapNodes(int node0, int node1) {
    std::swap(order_nodes[node0]->order, order_nodes[node1]->order);
    std::swap(order_nodes[node0], order_nodes[node1]);
}

void Graph::makeNodeInvisible(int order_of_node) {
    Node* node = order_nodes[order_of_node];
    for (int i = node->offset_visible_nodes; i < node->neighbours.size(); i++) {
        Node* neighbour = &fixed_nodes[node->neighbours[i]];

        //make node invisible in all adjacencie lists of fixed neighbours
        for (int j = neighbour->offset_visible_nodes; j < neighbour->neighbours.size(); j++) {
            if (neighbour->neighbours[j] != node->id) {
                continue;
            }
            //add node to invisible nodes and add visible counter
            std::swap(neighbour->neighbours[j], neighbour->neighbours[neighbour->offset_visible_nodes]);
            neighbour->offset_visible_nodes++;
        }
    }

    //set visible counter to end of adjacency array -> as if all edges removed
    node->offset_visible_nodes = node->neighbours.size();
}

std::pair<std::vector<Node*>, int> Graph::Greedy() {
    int crossings;
    for (int i = 0; i < order_nodes.size() - 1; i++) {
        int crossings_old = this->countCrossings();
        this->swapNodes(i, i + 1);
        int crossings_new = this->countCrossings();
        if (crossings_old < crossings_new) {
            this->swapNodes(i, i + 1);
        }
        if (i == order_nodes.size() - 2) {
            crossings = std::min(crossings_new, crossings_old);
        }
    }

    return std::make_pair(this->order_nodes, crossings);
}

void Graph::Median_Heuristic()
{
    for (int i = 0; i < order_nodes.size(); i++) {
        for (int j = 0; j < order_nodes[i]->neighbours.size(); j++) {
            order_nodes[i]->median += order_nodes[i]->neighbours[j];
        }

        order_nodes[i]->median = order_nodes[i]->median / order_nodes[i]->neighbours.size();
        std::cout << "id: " << order_nodes[i]->id << " median: " << order_nodes[i]->median << std::endl;
    }

    std::sort(order_nodes.begin(), order_nodes.end(), [](const Node* a, const Node* b) {
        return a->median < b->median;
        });

    for (int i = 0; i < order_nodes.size(); i++) {
        order_nodes[i]->order = i;
    }
}


bool Graph::verifier(Graph check)
{
    std::vector<bool> unique_id_vec(this->order_nodes.size(), false);

    for (int i = 0; i < order_nodes.size(); i++) {

        //check every node id comes exactly once
        if (unique_id_vec[this->order_nodes[i]->id - this->n0 - 1]) { return false; }
        unique_id_vec[this->order_nodes[i]->id - this->n0 - 1] = true;

        //check for "illegal node id's"
        if (this->order_nodes[i]->id > (check.n0 + check.n1)) { return false; }

        //neighbours of each node are the same as the graph at the beginning
        for (int j = 0; j < order_nodes[i]->neighbours.size(); j++) {
            int id_check = this->order_nodes[i]->id;
            if (this->order_nodes[i]->neighbours[j] != check.movable_nodes[id_check - this->n0 - 1].neighbours[j]) { return false; }
        }
    }

    //a node id disappeared from the order_nodes array.
    for (int i = 0; i < unique_id_vec.size(); i++) {
        if (!unique_id_vec[i]) { return false; }
    }
    return true;
}

std::pair<int, int> Graph::DFSforPartition(int start_node_fixed_id, int partition ,std::vector<bool>& visited) {    
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

        if (current_node < n0 && current_node > last_node_interval) {
            last_node_interval = current_node;
        }

        if (!visited[current_node]) {
            std::cout << "Visited " << current_node << " ";
            visited[current_node] = true;
            if (current_node >= n0) {
                // movable_nodes[current_node - n0].partition = partition;
                movable_nodes[current_node - n0].partition.push_back(partition);
            }
            else {
                fixed_nodes[current_node].partition.push_back(partition);
            }
        }

        if (current_node >= n0) {
            for (int i = 0; i < movable_nodes[current_node - n0].neighbours.size(); i++) {
                if (!visited[movable_nodes[current_node - n0].neighbours[i] - 1]) {
                    stack.push(movable_nodes[current_node - n0].neighbours[i] - 1);
                }
            }
        }
        else {
            for (int i = 0; i < fixed_nodes[current_node].neighbours.size(); i++) {
                if (!visited[fixed_nodes[current_node].neighbours[i] - 1]) {
                    stack.push(fixed_nodes[current_node].neighbours[i] - 1);
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

    std::vector<bool> visited(this->n0 + this->n1, false);

    while (rounds_counter < this->n0) {
        
        std::pair<int, int> end_start_node_interval = DFSforPartition(start_node_fix, partition, visited);
        rounds_counter++;
        int next_fixed_node = end_start_node_interval.first + 1;
        
        while (next_fixed_node <= end_start_node_interval.second) {
            if (visited[next_fixed_node]) {
                next_fixed_node++;
                rounds_counter++;
            }
            else {
                std::pair<int,int> end_start_node_interval_2 = DFSforPartition(next_fixed_node, partition, visited);
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

    if (start_node_id >= n0) {
        // Go through all vertices adjacent to this
        for (int i = 0; i < movable_nodes[start_node_id - n0].neighbours.size(); i++) {
            // If v is not visited yet, then make it a child of u
            // in DFS tree and recur for it
            if (!visited[movable_nodes[start_node_id - n0].neighbours[i] - 1]) {
                children++;
                APUtil(movable_nodes[start_node_id - n0].neighbours[i] - 1, visited, disc, low, time, start_node_id, isAP);

                // Check if the subtree rooted with v has
                // a connection to one of the ancestors of u
                low[start_node_id] = std::min(low[start_node_id], low[movable_nodes[start_node_id - n0].neighbours[i] - 1]);

                // If start_node_id is not root and low value of one of
                // its child is more than discovery value of u.
                if (parent != -1 && low[movable_nodes[start_node_id - n0].neighbours[i] - 1] >= disc[start_node_id])
                    isAP[start_node_id] = true;
            }

            // Update low value of start_node_id for parent function calls.
            else if (movable_nodes[start_node_id - n0].neighbours[i] - 1 != parent)
                low[start_node_id] = std::min(low[start_node_id], disc[movable_nodes[start_node_id - n0].neighbours[i] - 1]);
        }
    }
    else {
        // Go through all vertices adjacent to this
        for (int i = 0; i < fixed_nodes[start_node_id].neighbours.size(); i++) {
            // If v is not visited yet, then make it a child of u
            // in DFS tree and recur for it
            if (!visited[fixed_nodes[start_node_id].neighbours[i] - 1]) {
                children++;
                APUtil(fixed_nodes[start_node_id].neighbours[i] - 1, visited, disc, low, time, start_node_id, isAP);

                // Check if the subtree rooted with v has
                // a connection to one of the ancestors of u
                low[start_node_id] = std::min(low[start_node_id], low[fixed_nodes[start_node_id].neighbours[i] - 1]);

                // If start_node_id is not root and low value of one of
                // its child is more than discovery value of u.
                if (parent != -1 && low[fixed_nodes[start_node_id].neighbours[i] - 1] >= disc[start_node_id])
                    isAP[start_node_id] = true;
            }

            // Update low value of start_node_id for parent function calls.
            else if (fixed_nodes[start_node_id].neighbours[i] - 1 != parent)
                low[start_node_id] = std::min(low[start_node_id], disc[fixed_nodes[start_node_id].neighbours[i] - 1]);
        }
    }

    // If start_node_id is root of DFS tree and has two or more children.
    if (parent == -1 && children > 1)
        isAP[start_node_id] = true;
}

void Graph::AP()
{
    std::vector<int> disc(this->n0 + this->n1, 0);
    std::vector<int> low(this->n0 + this->n1);
    std::vector<int> visited(this->n0 + this->n1, false);
    std::vector<int> isAP(this->n0 + this->n1, false);
    int time = 0, par = -1;

    // Adding this loop so that the
    // code works even if we are given
    // disconnected graph
    for (int start_node_id= 0; start_node_id < this->n0 + this->n1; start_node_id++)
        if (!visited[start_node_id]) {
            APUtil(start_node_id, visited, disc, low, time, par, isAP);
        }

    // Printing the APs
    for (int start_node_id = 0; start_node_id < this->n0 + this->n1; start_node_id++) {
        if (isAP[start_node_id]) {
            std::cout << "Node id: " << start_node_id + 1 << " is AP" << std::endl;
        }
    }

    int partition_offset = 0;
    for (int i = 0; i < this->n0; i++) {
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
            if (g->getNodeByOrder(i)->neighbours.size() != g->getNodeByOrder(j)->neighbours.size()) {
                continue;
            }

            for (int k = 0; k < g->getNodeByOrder(i)->neighbours.size(); k++) {
                if (g->getNodeByOrder(i)->neighbours[k] != g->getNodeByOrder(j)->neighbours[k]) {
                    break;
                }

                //found twins
                if (k == g->getNodeByOrder(i)->neighbours.size()-1) {
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
    for (int a = 0; a < order_nodes.size()-1; a++) {
        for (int b = a+1; b < order_nodes.size(); b++) {
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
    std::vector<Node*> baseOrder = g.getOrderNodes();
    std::sort(baseOrder.begin(), baseOrder.end(), compareNodePointers);

    std::vector<Node*> bestOrder = baseOrder;
    int bestCrossings = g.countCrossings();

    do {
        Graph tmp = g;
        tmp.setOrderNodes(baseOrder);
        // Display the current permutation
        /*std::cout << "Permutation:";
        for (const Node* node : baseOrder) {
            std::cout << " order_nodes: " << node->order << ", ID: " << node->id;
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

int factorial(int n) {
    return (n == 0 || n == 1) ? 1 : n * factorial(n - 1);
}

std::pair<std::vector<Node*>, int> bruteForceParallel(Graph g) {
    std::vector<Node*> baseOrder = g.getOrderNodes();
    std::sort(baseOrder.begin(), baseOrder.end(), compareNodePointers);

    // Save current max
    std::vector<Node*> bestOrder = baseOrder;
    int bestCrossings = g.countCrossings();

    #pragma omp parallel for
    for (int i = 0; i < factorial(baseOrder.size()); i++){
        std::vector<Node*> permutation = baseOrder;

        // Generate the ith permutation
        for (int j = 0; j < i; ++j) {
            std::next_permutation(permutation.begin(), permutation.end(), compareNodePointers);
        }

        Graph tmp = g;
        tmp.setOrderNodes(permutation);

        int crossings = tmp.countCrossings();
        //std::cout << "Crossings: " << crossings << std::endl;
        #pragma omp critical
        {
            if (crossings < bestCrossings) {
                bestOrder = baseOrder;
                bestCrossings = crossings;
            }
        }
    }

    return std::make_pair(bestOrder, bestCrossings);
}