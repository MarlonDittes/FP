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
#include <queue>



//Constructor dependent on size of movable_nodes
Graph::Graph(int n0, int n1, int m) : n0(n0), n1(n1), m(m), graph(n0 + n1), order_nodes(n1), partitions(1) {

    for (int i = 0; i < n0; i++) {
        graph[i].id = i;
    }
    //initialize id and order in the graph
    for (int i = 0; i < n1; i++) {
        graph[n0 + i].id = n0 + i;
        graph[n0 + i].order = i;
        order_nodes[i] = &graph[n0 + i]; //link order to index (not id!)
    }
}

//adding Edges from movable_nodes to neighbours
void Graph::addEdge(int x, int y) {
    graph[x].neighbours.push_back(y);
    graph[y].neighbours.push_back(x);
}

void Graph::printGraph() {
    std::cout << "Printing graph..." << std::endl;
    std::cout << "fixed nodes = " << n0 << " moveable nodes = " << n1 << " total edges = " << m << std::endl;

    //print in order_nodes of movable_nodes (B)
    for (int i = 0; i < order_nodes.size(); i++) {
        std::cout << order_nodes[i]->id+1 << " -> ";
        for (int j = 0; j < order_nodes[i]->neighbours.size(); j++) {
            std::cout << order_nodes[i]->neighbours[j] +1<< " ";  
        }
        std::cout << std::endl;
    }
    /*
    for (int i = 0; i < graph.size(); i++) {
        std::cout << graph[i].id << " partition :  ";
        for (int j = 0; j < graph[i].partition.size(); j++) {
            std::cout << graph[i].partition[j] << " ";
        }
        std::cout << std::endl;
    }*/

    for (int i = 0; i < this->partitions.size(); i++) {
        std::cout << "partition : " << i <<" node in partition : ";
        for (int j = 0; j < this->partitions[i].size(); j++) {
            std::cout << this->partitions[i][j]->id +1<<" ";
        }
        std::cout << std::endl;
    }
}


void Graph::sortNeighbours() {
    for (int i = 0; i < graph.size(); i++) {
        std::sort(graph[i].neighbours.begin(), graph[i].neighbours.end());
    }
}

long Graph::countCrossingsMarlon() {
    if (m == 0) {
        return 0;
    }
    long crossings = 0;
    //iterate through visible movable_nodes
    for (int u = offset_visible_order_nodes; u < order_nodes.size() - 1; u++) {
        //iterate through visible movable_nodes+1
        for (int v = u + 1; v < order_nodes.size(); v++) {
            for (auto& u_neighbour : order_nodes[u]->neighbours) {
                for (auto& v_neighbour : order_nodes[v]->neighbours) {
                    if (u_neighbour > v_neighbour) {
                        crossings++;
                    }
                }
            }
        }
    }
    return crossings;
} 

long Graph::countCrossings() {
    if (m == 0) {
        return 0;
    }
    long crossings = 0;
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

long Graph::countCrossingsBranching() {
    if (m == 0) {
        return 0;
    }
    long crossings = 0;
    //iterate through movable_nodes
    for (int u = 0; u < order_nodes.size() - 1; u++) {
        //iterate through movable_nodes+1
        for (int v = u + 1; v < order_nodes.size(); v++) {
            for (int i = order_nodes[u]->offset_visible_nodes; i < order_nodes[u]->neighbours.size(); i++) {
                for (int j = order_nodes[v]->offset_visible_nodes; j < order_nodes[v]->neighbours.size(); j++) {
                    if (order_nodes[u]->neighbours[i] > order_nodes[v]->neighbours[j]) {
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

void Graph::swapNodesBranching(int node0, int node1) {
    std::swap(order_nodes[node0]->order, order_nodes[node1]->order);
    std::swap(order_nodes[node0], order_nodes[node1]);
    std::swap(graph[order_nodes[node0]->id].order, graph[order_nodes[node1]->id].order);
}

void Graph::makeNodeInvisibleMarlon(int order_of_node) {
    assert(offset_visible_order_nodes <= order_of_node);
    assert(order_of_node <= order_nodes.size() - 1);
    //Make sure node is not already invisible
    assert(order_nodes[order_of_node]->offset_visible_nodes != order_nodes[order_of_node]->neighbours.size());

    Node* node = order_nodes[order_of_node];
    for (int i = node->offset_visible_nodes; i < node->neighbours.size(); i++) {
        Node* neighbour = &graph[node->neighbours[i]];

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

    //move invisible node to beginning of order_nodes
    while(node->order != this->offset_visible_order_nodes){
        //swap with node one in front
        this->swapNodes(node->order, node->order - 1);
    }
    this->offset_visible_order_nodes++;
}

void Graph::makeNodeInvisible(int order_of_node) {
    /*Node* node = order_nodes[order_of_node];
    for (int i = node->offset_visible_nodes; i < node->neighbours.size(); i++) {
        Node* neighbour = &graph[node->neighbours[i]];

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
    node->offset_visible_nodes = node->neighbours.size();*/
}

void Graph::makeNodeInvisibleBranching(int order_of_node) {
    Node* node = order_nodes[order_of_node];
    for (int i = node->offset_visible_nodes; i < node->neighbours.size(); i++) {
        Node* neighbour = &graph[node->neighbours[i]];

        //make node invisible in all adjacencie lists of fixed neighbours
        for (int j = neighbour->offset_visible_nodes; j < neighbour->neighbours.size(); j++) {
            if (neighbour->neighbours[j] != node->id) {
                continue;
            }
            //add node to invisible nodes and asdd visible counter
            std::swap(neighbour->neighbours[j], neighbour->neighbours[neighbour->offset_visible_nodes]);
            neighbour->offset_visible_nodes++;
        }
    }

    //set visible counter to end of adjacency array -> as if all edges removed
    node->offset_visible_nodes = node->neighbours.size();
    graph[node->id].offset_visible_nodes = node->neighbours.size();
}

void Graph::makeNodeVisibleMarlon() {
    assert(offset_visible_order_nodes >= 1);

    int order_of_node = this->offset_visible_order_nodes - 1;
    //Make sure node is really invisible
    assert(order_nodes[order_of_node]->offset_visible_nodes == order_nodes[order_of_node]->neighbours.size());

    Node* node = order_nodes[order_of_node];
    node->offset_visible_nodes = 0;

    for (int i = 0; i < node->neighbours.size(); i++) {
        Node* neighbour = &graph[node->neighbours[i]];

        //make node visible in all adjacencie lists of fixed neighbours
        for (int j = 0; j < neighbour->offset_visible_nodes; j++) {
            if (neighbour->neighbours[j] != node->id) {
                continue;
            }
            //add node to visible nodes and reduce visible counter
            std::swap(neighbour->neighbours[j], neighbour->neighbours[neighbour->offset_visible_nodes - 1]);
            neighbour->offset_visible_nodes--;
        }
    }

    this->offset_visible_order_nodes--;
}

void Graph::makeNodeVisible(int order_of_node) {
    Node* node = order_nodes[order_of_node];
    node->offset_visible_nodes = 0;
    graph[node->id].offset_visible_nodes = 0;

    for (int i = 0; i < node->neighbours.size(); i++) {
        Node* neighbour = &graph[node->neighbours[i]];

        //make node visible in all adjacencie lists of fixed neighbours
        for (int j = 0; j < neighbour->offset_visible_nodes; j++) {
            if (neighbour->neighbours[j] != node->id) {
                continue;
            }
            //add node to visible nodes and reduce visible counter
            std::swap(neighbour->neighbours[j], neighbour->neighbours[neighbour->offset_visible_nodes - 1]);
            neighbour->offset_visible_nodes--;
        }
    }
}

std::pair<std::vector<Node*>, long> Graph::Greedy() {
    long crossings;
    for (int i = 0; i < order_nodes.size() - 1; i++) {
        int crossings_old = this->countCrossings();
        this->swapNodes(i, i + 1);
        int crossings_new = this->countCrossings();
        if (crossings_old < crossings_new) {
            this->swapNodes(i, i + 1);
        }

        // Case for last comparison
        if (i == order_nodes.size() - 2) {
            crossings = std::min(crossings_new, crossings_old);
        }
    }

    return std::make_pair(this->order_nodes, crossings);
}

void Graph::MedianHeuristicMarlon() {
    for (int i = offset_visible_order_nodes; i < order_nodes.size(); i++) {
        order_nodes[i]->median = 0;
        for (int j = order_nodes[i]->offset_visible_nodes; j < order_nodes[i]->neighbours.size(); j++) {
            order_nodes[i]->median += order_nodes[i]->neighbours[j];
        }
        if ((order_nodes[i]->neighbours.size() - order_nodes[i]->offset_visible_nodes) != 0){
            order_nodes[i]->median = order_nodes[i]->median / (order_nodes[i]->neighbours.size() - order_nodes[i]->offset_visible_nodes);
            //std::cout << "id: " << order_nodes[i]->id << " median: " << order_nodes[i]->median << std::endl;
        } else {
            order_nodes[i]->median = 0;
        }
        
    }

    std::sort(order_nodes.begin() + offset_visible_order_nodes, order_nodes.end(), [](const Node* a, const Node* b) {
        return a->median < b->median;
        });

    for (int i = 0; i < order_nodes.size(); i++) {
        order_nodes[i]->order = i;
    }
}

void Graph::Median_Heuristic() {
    for (int i = 0; i < order_nodes.size(); i++) {
        for (int j = 0; j < order_nodes[i]->neighbours.size(); j++) {
            order_nodes[i]->median += order_nodes[i]->neighbours[j];
        }

        order_nodes[i]->median = order_nodes[i]->median / order_nodes[i]->neighbours.size();
        //std::cout << "id: " << order_nodes[i]->id << " median: " << order_nodes[i]->median << std::endl;
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
        if (unique_id_vec[this->order_nodes[i]->id - this->n0]) { return false; }
        unique_id_vec[this->order_nodes[i]->id - this->n0] = true;

        //check for "illegal node id's"
        if (this->order_nodes[i]->id >= (check.n0 + check.n1)) { return false; }

        //neighbours of each node are the same as the graph at the beginning
        for (int j = 0; j < order_nodes[i]->neighbours.size(); j++) {
            int id_check = this->order_nodes[i]->id;
            if (this->order_nodes[i]->neighbours[j] != check.graph[id_check].neighbours[j]) { return false; }
        }
    }

    //a node id disappeared from the order_nodes array.
    for (int i = 0; i < unique_id_vec.size(); i++) {
        if (!unique_id_vec[i]) { return false; }
    }
    return true;
}

bool compareNodeID(Node* a, Node* b) {
    return a->id < b->id;
}

bool compareNodeOrder(Node* a, Node* b) {
    return a->order < b->order;
}

//change function such that it works for degree 0 nodes.
std::pair<int, int> Graph::DFSforPartition(int start_node_fixed_id, int partition, std::vector<bool>& visited) {
    // need to make sure at the beginning that all nodes with degree 0 are removed. we can do this since by setting their
    // order to be the right most position and adding an ignore marker.
    Stack stack;

    stack.push(start_node_fixed_id);

    bool already_in_partition = false;
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
            //std::cout << "Visited " << current_node << " ";
            visited[current_node] = true;
            this->graph[current_node].partition.push_back(partition);
            for (int i = 0; i < graph[current_node].partition.size(); i++) {
                if (partition == graph[current_node].partition[i]) {
                    already_in_partition = true;
                }
            }
            if (!already_in_partition) {
                this->graph[current_node].partition.push_back(partition);
                already_in_partition = false;
            }
        }

        for (int i = graph[current_node].offset_visible_nodes; i < graph[current_node].neighbours.size(); i++) {
            if (!visited[graph[current_node].neighbours[i]]) {
                stack.push(graph[current_node].neighbours[i]);
            }
        }
    }

    return std::make_pair(start_node_interval, last_node_interval);
}

void Graph::Partition()
{
    int start_node_fix = graph[0].id;
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
                std::pair<int, int> end_start_node_interval_2 = DFSforPartition(next_fixed_node, partition, visited);
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
void Graph::APUtil(int start_node_id, std::vector<bool>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<int>& isAP)
{
    // Count of children in DFS Tree
    int children = 0;

    // Mark the current node as visited
    visited[start_node_id] = true;

    // Initialize discovery time and low value
    disc[start_node_id] = low[start_node_id] = ++time;

    // Go through all vertices adjacent to this
    for (int i = graph[start_node_id].offset_visible_nodes; i < graph[start_node_id].neighbours.size(); i++) {
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (!visited[graph[start_node_id].neighbours[i]]) {
            children++;
            APUtil(graph[start_node_id].neighbours[i], visited, disc, low, time, start_node_id, isAP);

            // Check if the subtree rooted with v has
            // a connection to one of the ancestors of u
            low[start_node_id] = std::min(low[start_node_id], low[graph[start_node_id].neighbours[i]]);

            // If start_node_id is not root and low value of one of
            // its child is more than discovery value of u.
            if (parent != -1 && low[graph[start_node_id].neighbours[i]] >= disc[start_node_id])
                isAP[start_node_id] = true;
        }

        // Update low value of start_node_id for parent function calls.
        else if (graph[start_node_id].neighbours[i] != parent)
            low[start_node_id] = std::min(low[start_node_id], disc[graph[start_node_id].neighbours[i]]);
    }

    // If start_node_id is root of DFS tree and has two or more children.
    if (parent == -1 && children > 1)
        isAP[start_node_id] = true;
}

//change function such that it works for degree 0 nodes.
void Graph::CreatePartitionsVector(int start_node_fixed_id, int& partition_id, std::vector<bool>& visited) {
    // need to make sure at the beginning that all nodes with degree 0 are removed. we can do this since by setting their
    // order to be the right most position and adding an ignore marker.
    Stack stack;
    stack.push(start_node_fixed_id);

    if (partition_id >= this->partitions.size()) {
        std::vector<Node*> partition;
        this->partitions.push_back(partition);
    }

    bool partition_id_in_node = false;
    bool partition_in_partition_vector = false;

    //move this to a seperate function called DFS
    while (!stack.isEmpty()) {
        int current_node = stack.peek();
        stack.pop();

        if (!visited[current_node]) {
            visited[current_node] = true;

            if (graph[current_node].partition.size() > 1) {
                partition_id++;
            }

            //change the partition vector in the node data structure if necessary
            for (int i = 0; i < graph[current_node].partition.size(); i++) {
                if (graph[current_node].partition[i] == partition_id) {
                    partition_id_in_node = true;
                }
            }

            if (!partition_id_in_node) {
                if (graph[current_node].isAP) {
                    graph[current_node].partition.push_back(partition_id);
                }
                else {
                    graph[current_node].partition[graph[current_node].partition.size() - 1] = partition_id;
                }
            }

            //add node to partitionsvector
            for (int i = 0; i < partitions[partition_id].size(); i++) {
                if (partitions[partition_id][i]->id == current_node) {
                    partition_in_partition_vector = true;
                }
            }
            if (!partition_in_partition_vector) {
                this->partitions[partition_id].push_back(&graph[current_node]);
            }
        }

        //add neighbours
        for (int i = graph[current_node].offset_visible_nodes; i < graph[current_node].neighbours.size(); i++) {
            if (!visited[graph[current_node].neighbours[i]]) {
                    if(graph[graph[current_node].neighbours[i]].isAP) {
                        //just add the AP point to the partitions vector but not to the working stack
                        partition_in_partition_vector = false;
                        for (int j = 0; j < partitions[partition_id].size(); j++) {
                            if (partitions[partition_id][j]->id == graph[current_node].neighbours[i]) {
                                partition_in_partition_vector = true;
                            }
                        }
                        if (!partition_in_partition_vector) {
                            this->partitions[partition_id].push_back(&graph[graph[current_node].neighbours[i]]);
                        }

                        int node_partition_size = graph[graph[current_node].neighbours[i]].partition.size();
                        graph[graph[current_node].neighbours[i]].partition[node_partition_size - 1] = partition_id;
                    }
                    else {
                        stack.push(graph[current_node].neighbours[i]);
                    }
            }
        }
    }
}


void Graph::AP()
{
    std::vector<int> disc(this->n0 + this->n1, 0);
    std::vector<int> low(this->n0 + this->n1);
    std::vector<bool> visited(this->n0 + this->n1, false);
    std::vector<int> isAP(this->n0 + this->n1, false);
    int time = 0, par = -1;

    // Adding this loop so that the
    // code works even if we are given
    // disconnected graph
    for (int start_node_id = 0; start_node_id < this->n0 + this->n1; start_node_id++)
        if (!visited[start_node_id]) {
            APUtil(start_node_id, visited, disc, low, time, par, isAP);
        }

    visited = std::vector<bool>(this->n0 + this->n1, false);
    std::vector<bool> visited_node_stack(this->n0 + this->n1, false);
    int partition_id = 0;

    for (int i = 0; i < n0; i++) {
        if (isAP[i]) { 
            std::cout << "node : " << i << " is AP" << std::endl;
            graph[i].isAP = true;
        }
    }

    Stack stack;
    stack.push(n0);

    //move this to a seperate function called DFS
    while (!stack.isEmpty()) {
        int current_node = stack.peek();
        stack.pop();

        if (!visited[current_node]) {
            CreatePartitionsVector(current_node, partition_id, visited);
            partition_id++;
        }

        for (int i = graph[current_node].offset_visible_nodes; i < graph[current_node].neighbours.size(); i++) {
            if (!visited[graph[current_node].neighbours[i]]) {
                stack.push(graph[current_node].neighbours[i]);
            }
        }

        if (stack.isEmpty()) {
            for (int i = 0; i < n0 + n1; i++) {
                if (!visited[graph[i].id] && graph[i].offset_visible_nodes != graph[i].neighbours.size()) {
                    stack.push(i);
                    break;
                }
            }
        }
    }

    /*CreatePartitionsVector(n0, partition_id, visited);
    partition_id++;
    CreatePartitionsVector(1, partition_id, visited);
    partition_id++;

    for (int start_node_id = 0; start_node_id < this->n0 + n1; start_node_id++) {
        if (!visited[start_node_id]) {
            CreatePartitionsVector(start_node_id, partition_id, visited);
            partition_id++;
        }
    }*/

    for (int i = 0; i < partitions.size(); i++) {
        std::sort(partitions[i].begin(), partitions[i].end(), compareNodeID);
    }

    for (int i = 0; i < this->n0 + this->n1; i++) {
        graph[i].isAP = false;
    }
}


bool Graph::DFS_for_sorted_straight_line(int start_node, std::vector<bool>& visited) {
    // need to make sure at the beginning that all nodes with degree 0 are removed. we can do this since by setting their
    // order to be the right most position and adding an ignore marker.
    Stack stack;
    stack.push(start_node);
    int mod_calc = 0;
    int last_fix_node_id = start_node;


    //move this to a seperate function called DFS
    while (!stack.isEmpty()) {

        int current_node = stack.peek();
        stack.pop();

        if (mod_calc % 2 == 0 && current_node != 0) {
            if (last_fix_node_id != current_node - 1) {
                std::cout << "The given Graph / Partition is not a straight line" << std::endl;
                return false;
            }
            else {
                last_fix_node_id = current_node;
            }
        }

        if (!visited[current_node]) {
            std::cout << "Visited " << current_node << " ";
            visited[current_node] = true;
        }
        for (int i = 0; i < graph[current_node].neighbours.size(); i++) {
            if ((graph[current_node].neighbours.size() - graph[current_node].offset_visible_nodes) > 2) {
                std::cout << "node : " << current_node << " has more than 1 neighbour therefore straight line reduction can not be achieved " << std::endl;
                return false;
            }
            if (!visited[graph[current_node].neighbours[i]]) {
                stack.push(graph[current_node].neighbours[i]);
            }
        }

        mod_calc++;
    }

    return true;
}



void Graph::Sorted_straight_line_reduction() {
    // max 1 DFS
    // after every second search make sure the node that is visited has an ID exactly +1 of the node before.
    // this can be done using a mod 2 if statement to check if the node we are about the visit is a fix node or not
    // all nodes must have a max of 1 neighbour. to make sure of this after 1 DFS check all nodes were visited.

    // after we have found out that it is a straight line input that can lead to 0 crossings --> use median Heuristic,
    // to calculate the median of their neighbours, which then in turn gives us the order of the moveable nodes.
    int start_node_id = 0;
    std::vector<bool> visited(n0 + n1);
    if (DFS_for_sorted_straight_line(start_node_id, visited)) {
        for (int i = 0; i < visited.size(); i++) {
            if (!visited[i]) {
                std::cout << "not a sorted straight line reduction, not all nodes in this given partition/graph can be visited" << std::endl;
            }
        }
        Median_Heuristic();
    }
    else {
        std::cout << "sorted straight line reduction could not be done" << std::endl;
    }
}

void Graph::setNoPartitioning() {
    if (this->partitions.empty()){
        std::vector<Node*> tmp(0);
        for (auto& node : this->graph){
            tmp.push_back(&node);
        }
        this->partitions.push_back(tmp);
    }
}

void Graph::setOrderNodes(std::vector<Node*> order){
    this->order_nodes = order;

    for (int i = offset_visible_order_nodes; i < order_nodes.size(); i++){
        order_nodes[i]->order = i;
    }
}

void Branch_and_Bound(Graph* G) {
    G->Median_Heuristic();
    Graph verifier = *G;
    std::vector<Node*>best_configuration = G->getOrderNodes();
    int best_solution = G->countCrossingsBranching();

    std::cout << "Median Heuristic : " << best_solution << std::endl;

    Graph G_branch = *G;
    for (int i = 1; i < G_branch.getOrderNodes().size(); i++) {
        G_branch.makeNodeInvisibleBranching(i);
    }
    
    exploreBranch(G_branch, *G, 0, best_solution, best_configuration);

    std::cout << "best solution "<<best_solution << std::endl;
    for (int i = 0; i < best_configuration.size(); i++) {
        std::cout << best_configuration[i]->id << std::endl;
    }

    G->setOrderNodes(best_configuration);
    std::cout << "num of crossings : " << G->countCrossings() << std::endl;

    if (G->verifier(verifier)) {
        std::cout << "Graph is valid" << std::endl;
    }
    else {
        std::cout << "Graph is NOT valid" << std::endl;
    }
}

void exploreBranch(Graph G_branch, Graph& G_original, int depth, int& best_solution, std::vector<Node*>& best_configuration) {

    if (depth == G_original.getN1() - 1) {
        int current_solution = G_branch.countCrossingsBranching();
        if (current_solution < best_solution) {
            best_solution = current_solution;
            best_configuration = G_branch.getOrderNodes();
        }
        return;
    }

    G_branch.makeNodeVisible(depth + 1);

    if (depth == 0) {
        if (G_branch.countCrossingsBranching() <= best_solution) {
            exploreBranch(G_branch, G_original, depth + 1, best_solution, best_configuration);
        }
        G_branch.swapNodesBranching(0, 1);
        if (G_branch.countCrossingsBranching() <= best_solution) {
            exploreBranch(G_branch, G_original, depth + 1, best_solution, best_configuration);
        }
    }
    else {
        if (G_branch.countCrossingsBranching() <= best_solution) {
            exploreBranch(G_branch, G_original, depth + 1, best_solution, best_configuration);
        }
        for (int i = depth + 1; i > 0; i--) {
            G_branch.swapNodesBranching(i, i - 1);
            if (G_branch.countCrossingsBranching() <= best_solution) {
                exploreBranch(G_branch, G_original, depth + 1, best_solution, best_configuration);
            }
        }
    }
}



typedef std::vector<std::pair<Node*, Node*>> Twins;

Twins findTwins(Graph* g) {
    Twins twins;
    for (int i = 0; i < g->getSizeOfOrder() - 1; i++) {
        for (int j = i + 1; j < g->getSizeOfOrder(); j++) {
            if (g->getNodeByOrder(i)->neighbours.size() != g->getNodeByOrder(j)->neighbours.size()) {
                continue;
            }

            for (int k = 0; k < g->getNodeByOrder(i)->neighbours.size(); k++) {
                if (g->getNodeByOrder(i)->neighbours[k] != g->getNodeByOrder(j)->neighbours[k]) {
                    break;
                }

                //found twins
                if (k == g->getNodeByOrder(i)->neighbours.size() - 1) {
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

std::pair<std::vector<Node*>, long> bruteForce(Graph* g) {
    std::vector<Node*> baseOrder = g->getOrderNodes();
    std::sort(baseOrder.begin(), baseOrder.end(), compareNodeID);

    std::vector<Node*> bestOrder = baseOrder;
    long bestCrossings = g->countCrossings();

    do {
        g->setOrderNodes(baseOrder);
        // Display the current permutation
        /*std::cout << "Permutation:";
        for (const Node* node : baseOrder) {
            std::cout << " order_nodes: " << node->order << ", ID: " << node->id;
        }
        std::cout << "\n";
*/
        long crossings = g->countCrossings();
        //std::cout << "Crossings: " << crossings << std::endl;

        if (crossings < bestCrossings) {
            bestOrder = baseOrder;
            bestCrossings = crossings;
        }

    } while (std::next_permutation(baseOrder.begin(), baseOrder.end(), compareNodeID));

    return std::make_pair(bestOrder, bestCrossings);
}

std::pair<std::vector<Node*>, long> bruteForceOnSubgraph(Graph* g, int begin, int end) {
    std::vector<Node*> baseOrder = g->getOrderNodes();
    std::sort(baseOrder.begin() + begin, baseOrder.begin() + end, compareNodeID);

    std::vector<Node*> bestOrder = baseOrder;
    long bestCrossings = g->countCrossings();

    do {
        g->setOrderNodes(baseOrder);
        // Display the current permutation
        /*std::cout << "Permutation:";
        for (const Node* node : baseOrder) {
            std::cout << " order_nodes: " << node->order << ", ID: " << node->id;
        }
        std::cout << "\n";
*/
        long crossings = g->countCrossings();
        //std::cout << "Crossings: " << crossings << std::endl;

        if (crossings < bestCrossings) {
            bestOrder = baseOrder;
            bestCrossings = crossings;
        }

    } while (std::next_permutation(baseOrder.begin() + begin, baseOrder.begin() + end, compareNodeID));

    return std::make_pair(bestOrder, bestCrossings);
}

Graph* createGraphByPartition(Graph* g, std::vector<Node*> partition) {
    std::sort(partition.begin(), partition.end(), compareNodeID);
    int n0 = 0;
    int n1 = 0;
    int m = 0;
    std::vector<int> partitionIndices(0);

    //Count n0, n1, m
    int originalN0 = g->getN0();
    for (auto& node : partition) {
        partitionIndices.push_back(node->id);
        
        if (node->id < originalN0) {
            n0++;
        }
        //Count edges only once on moveable nodes, since AP are fix
        else {
            n1++;
            for (int i = node->offset_visible_nodes; i < node->neighbours.size(); i++){
                m++;
            }
        }
    }
    
    Graph* partGraph = new Graph(n0, n1, m);
    
    for (int i = 0; i < partition.size(); i++) {
        partGraph->setOldID(i, partition[i]->id);
        //Need this so to not add forward AND backward edge
        if (partition[i]->id >= originalN0) {
            for (int j = partition[i]->offset_visible_nodes; j < partition[i]->neighbours.size(); j++){
                partGraph->addEdge(i, std::distance(partitionIndices.begin(), std::find(partitionIndices.begin(), partitionIndices.end(), partition[i]->neighbours[j])));
            }
        }
    }
    partGraph->sortNeighbours();
    return partGraph;
}

std::pair<std::vector<Node*>, long> branching (Graph* g, std::vector<Reduction> reductionTypes){
    bool changed = false;
    //for (auto reduct : reductionTypes) {
        /*
        Need to implement this somehow

        if (reduct.reduce()){
            reduct.apply();
            changed = true;
        }
        */
    //}
    
    //Reduce our instance if no more reductions applicable
    std::pair<std::vector<Node*>, long> result;
    if (!changed) {
        auto orderNodes = g->getOrderNodes();
        int visibleNodeOffset = g->getOffsetVisibleOrderNodes();

        if ((orderNodes.size() - visibleNodeOffset) > 2){
            // Find highest degree visible moveable node
            Node* maxDegreeNode = orderNodes[visibleNodeOffset];
            int maxDegree = orderNodes[visibleNodeOffset]->neighbours.size() - orderNodes[visibleNodeOffset]->offset_visible_nodes;
            for (int i = visibleNodeOffset + 1; i < orderNodes.size(); i++){
                if (maxDegree < orderNodes[i]->neighbours.size() - orderNodes[i]->offset_visible_nodes){
                    maxDegree = orderNodes[i]->neighbours.size() - orderNodes[i]->offset_visible_nodes;
                    maxDegreeNode = orderNodes[i];
                }
            }

            // Remove node for now
            g->makeNodeInvisibleMarlon(maxDegreeNode->order);
            // Solve on remaining nodes
            result = BranchAndReduce(g, reductionTypes);
            // Add node back
            g->makeNodeVisibleMarlon();

            // Calculate Crossings for every possible position of node
            int minCrossings = g->countCrossingsMarlon();
            auto minOrder = g->getOrderNodes();
            for (int i = visibleNodeOffset + 1; i < orderNodes.size(); i++){
                g->swapNodes(maxDegreeNode->order, i);
                int crossings = g->countCrossingsMarlon();
                if (crossings < minCrossings){
                    minCrossings = crossings;
                    minOrder = g->getOrderNodes();
                }
            }

            // Reconstruct min order
            g->setOrderNodes(minOrder);

            result = std::make_pair(g->getOrderNodes(), minCrossings);

        } else if ((orderNodes.size() - visibleNodeOffset) == 2){
            int crossings1 = g->countCrossingsMarlon();
            g->swapNodes(visibleNodeOffset, visibleNodeOffset + 1);
            int crossings2 = g->countCrossingsMarlon();
            // go back to previous order if it was better
            if (crossings1 < crossings2){
                g->swapNodes(visibleNodeOffset, visibleNodeOffset + 1);
                result = std::make_pair(g->getOrderNodes(), crossings1);
            } else {
                result = std::make_pair(g->getOrderNodes(), crossings2);
            }
 
        } else {
            // Less than two moveable nodes left -> nothing left to do
            //std::cout << "Only one moveable node."<< std::endl;
            result = std::make_pair(g->getOrderNodes(), 0);
        }      
    }
    else {
        result = BranchAndReduce(g, reductionTypes);
    }

    return result;
}

std::pair<std::vector<Node*>, long> BranchAndReduce(Graph* g, std::vector<Reduction> reductionTypes) {
    g->MedianHeuristicMarlon();

    g->Partition();
    g->AP();
    //g->setNoPartitioning();

    std::vector<std::vector<Node*>>& partitions = g->getPartitions();
    std::vector<Node*> solution(0);
    long sumCrossings = 0;
    if (partitions.size() > 1){
        std::vector<std::pair<std::vector<Node*>, long>> results(0);

        // get sub solutions
        for (auto& part : partitions) {
            Graph* partGraph = createGraphByPartition(g, part);

            auto result = branching(partGraph, reductionTypes);
            results.push_back(result);
        }

        // Combine sub solutions
        std::vector<Node*> subSolutions(0);
        for (auto& result : results) {
            subSolutions.insert(subSolutions.end(), result.first.begin(), result.first.end());
            sumCrossings += result.second;
        }

        // Apply solution to original graph
        std::vector<Node*> oldOrder = g->getOrderNodes();
        auto& nodes = g->getGraph();
        for (int i = 0; i < g->getOffsetVisibleOrderNodes(); i++){
            solution.push_back(oldOrder[i]);
        }
        for (auto& node : subSolutions){
            solution.push_back(&nodes[node->old_id]);
        }

        g->setOrderNodes(solution);

    } else {
        auto result = branching(g, reductionTypes);
        solution = result.first;
        sumCrossings = result.second;
    }
    
    return std::make_pair(solution, sumCrossings);
}