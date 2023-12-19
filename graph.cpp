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
Graph::Graph(int n0, int n1, int m) : n0(n0), n1(n1), m(m), graph(n0+n1), order_nodes(n1) {

    for (int i = 0; i < n0; i++) {
        graph[i].id = i;
    }
    //initialize id and order in the graph
    for (int i = 0; i < order_nodes.size(); i++) {
        graph[n0+i].id = n0+i; 
        graph[n0+i].order = i;
        order_nodes[i] = &graph[n0+i]; //link order to index (not id!)
    }
}

//adding Edges from movable_nodes to neighbours
void Graph::addEdge(int x, int y) {
    graph[x].neighbours.push_back(y);
    graph[y].neighbours.push_back(x);   
}

void Graph::printGraph() {
    std::cout << "Printing graph..." << std::endl;

    //print in order_nodes of movable_nodes (B)
    for (int i = 0; i < order_nodes.size(); i++) {
        std::cout << order_nodes[i]->id << " -> ";
        for (int j = 0; j < order_nodes[i]->neighbours.size(); j++) {
            std::cout << order_nodes[i]->neighbours[j] << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < graph.size(); i++) {
        std::cout << "Node Id: " << graph[i].id << " Partition : ";
        for (int j = 0; j < graph[i].partition.size(); j++) {
            std::cout << graph[i].partition[j] << " ";
        }
        std::cout << std::endl;
    }
}


void Graph::sortNeighbours() {
    for (int i = 0; i < graph.size(); i++) {
        std::sort(graph[i].neighbours.begin(), graph[i].neighbours.end());
    }
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

void Graph::Median_Heuristic()
{
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
            //std::cout << "Visited " << current_node << " ";
            visited[current_node] = true;
            this->graph[current_node].partition.push_back(partition);
            this->partitions[partition].push_back(current_node);
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

        std::vector<Node*> tmp(0);
        this->partitions.push_back(tmp);
        
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
    for (int start_node_id = 0; start_node_id < this->n0 + this->n1; start_node_id++)
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
    /*for (int i = 0; i < this->n0; i++) {
        for (int j = 0; j < fixed_nodes[i].partition[j]; j++) {
            fixed_nodes[i].partition[j] = fixed_nodes[i].partition[j] + partition_offset;

        }
        if (isAP[i]) {
            partition_offset++;


        }
    }*/

}


bool Graph::DFS_for_sorted_straight_line(int start_node,  std::vector<bool>& visited) {
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

        if(mod_calc % 2 == 0 && current_node != 0) { 
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
                return false; }
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


typedef std::vector<std::pair<Node*, Node*>> Twins;

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
    return a->id < b->id;
}

std::pair<std::vector<Node*>, long> bruteForce(Graph* g) {
    std::vector<Node*> baseOrder = g->getOrderNodes();
    std::sort(baseOrder.begin(), baseOrder.end(), compareNodePointers);

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

    } while (std::next_permutation(baseOrder.begin(), baseOrder.end(), compareNodePointers));

    return std::make_pair(bestOrder, bestCrossings);
}

std::pair<std::vector<Node*>, long> bruteForceOnSubgraph(Graph* g, int begin, int end){
    std::vector<Node*> baseOrder = g->getOrderNodes();
    std::sort(baseOrder.begin() + begin, baseOrder.begin() + end, compareNodePointers);

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

    } while (std::next_permutation(baseOrder.begin() + begin, baseOrder.begin() + end, compareNodePointers));

    return std::make_pair(bestOrder, bestCrossings);
}

// CHANGE THIS TO BE SOMEWHERE ELSE; THIS SHOULDNT BE HERE
enum Reduction {ZeroEdge, Complete};
constexpr int BRUTE_CUTOFF;

std::pair<std::vector<Node*>, long> BranchAndReduce (Graph* g, std::vector<Reduction> reductionTypes){
    g->Median_Heuristic();

    g->Partition();
    g->AP();

    auto partitions = g->getPartitions();           // NEED TO ADJUST THIS PROBABLY TO BE VECTOR<NODE*> NOT VECTOR<INT>!!, wait for shai AP
    std::vector<std::pair<std::vector<Node*>, long>> results(0);
    for (auto part : partitions){
        Graph* partGraph = createGraphByPartition(g, part);
        
        bool changed = false;
        for (auto reduct : reductionTypes){
            /*
            Need to implement this somehow

            if (reduct.reduce()){
                reduct.apply();
                changed = true;
            }
            */
        }

        //Check if need to brute force since no more reductions were applicable
        std::pair<std::vector<Node*>, long> result;

        //decide if we should do this:
        int numberOfVertices = g->getN0() + g->getN1();
        bool otherCondition = numberOfVertices <= BRUTE_CUTOFF;

        if (changed == false || otherCondition){
            result = bruteForce(g);
        } else {
            result = BranchAndReduce(partGraph, reductionTypes);
        }

        results.push_back(result);
    }

     std::vector<Node*> solution(0);
     long sumCrossings = 0;
     for (auto result : results){
        solution.insert(solution.end(), result.first.begin(), result.first.end());
        sumCrossings += result.second;
     }

     return std::make_pair(solution, sumCrossings);
}

Graph* createGraphByPartition(Graph* g, std::vector<int> partition){
    int n0 = 0;
    int n1 = 0;
    int m = 0;
    std::vector<Node> g_nodes = g->getGraph();
    std::vector<std::pair<int,int>> edges;
    for (auto& ind : partition){
        if (ind < g->getN0()){
            n0++;
            for (auto& neighbour : g_nodes[ind].neighbours){
                m++;
                edges.push_back(std::make_pair(ind, neighbour));
            }
        } else {
            n1++;
        }
    }

    Graph* partGraph = new Graph(n0,n1,m);
    for (auto& edge : edges){
        partGraph->addEdge(edge.first, edge.second);
    }

    return partGraph;
}

bool reduceCompleteReduction(Graph* g){
    int n0 = g->getN0();
    int n1 = g->getN1();
    int m = g->getM();

    if (n0 * n1 == m){
        return true;
    }
}

bool reduceZeroEdgeReduction(Graph* g){
    if (g->getM() == 0){
        return true;
    }
}