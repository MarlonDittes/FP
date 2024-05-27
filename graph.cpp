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
#include <stack>
#include <cstdlib> // For rand() and srand()
#include <ctime>   // For time()
#include <cstdint>
#include <string>

/*
#include "src_henning/src/definitions.h"
#include "src_henning/src/macros.h"
#include "src_henning/src/misc.h"
#include "src_henning/src/graph_hen.h"
#include "src_henning/src/solver_bf.h"
#include "src_henning/src/exhaustive_solver.h"
#include "src_henning/src/partitioner.h"
#include "src_henning/src/solver.h"
#include "src_henning/src/useless_reducer.h"
#include "src_henning/src/front_back_reducer.h"
*/

#include "TomAlv/heuristic_algorithm.h"
#include "TomAlv/heuristic_graph.h"
#include "TomAlv/median_algorithm.h"

//Constructor dependent on size of movable_nodes
Graph::Graph(int n0, int n1, int m) : n0(n0), n1(n1), m(m), activeEdges(m), graph(n0 + n1), order_nodes(n1), partitions(1) {

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
    graph[x].edges.push_back({ y, 1 });
    graph[y].edges.push_back({ x, 1 });
    graph[y].hash += x + 1; //add 1 because id starts at 0
}

void Graph::printGraph() {
    std::cout << "Printing graph..." << std::endl;
    std::cout << "fixed nodes = " << n0 << " moveable nodes = " << n1 << " total edges = " << m  << " offset_visible = " << offset_visible_order_nodes << std::endl;

    //print in order_nodes of movable_nodes (B)
    for (int i = 0; i < order_nodes.size(); i++) {
        std::cout << order_nodes[i]->id + 1 << " -> ";
        for (int j = 0; j < order_nodes[i]->edges.size(); j++) {
            std::cout << order_nodes[i]->edges[j].neighbour_id + 1 << " ";
        }
        std::cout << std::endl;
    }
    /*
    for (int i = 0; i < graph.size(); i++) {
        std::cout << graph[i].id + 1 << " partition :  ";
        for (int j = 0; j < graph[i].partition.size(); j++) {
            std::cout << graph[i].partition[j] << " ";
        }
        std::cout << std::endl;
    }
    */


    for (int i = 0; i < this->partitions.size(); i++) {
        std::cout << "partition : " << i << " node in partition : ";
        for (int j = 0; j < this->partitions[i].size(); j++) {
            std::cout << this->partitions[i][j]->id + 1 << " ";
        }
        std::cout << std::endl;
    }
}

void Graph::printGraphByPartitions() {
    for (int i = 0; i < partitions.size(); i++) {
        std::cout << "Partition " << i << ": " << std::endl;
        for (int j = 0; j < partitions[i].size(); j++) {
            std::cout << partitions[i][j]->id << " -> ";
            for (int k = 0; k < partitions[i][j]->edges.size(); k++) {
                std::cout << partitions[i][j]->edges[k].neighbour_id << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void Graph::sortOrderNodesByOrder() {
    std::sort(order_nodes.begin() + offset_visible_order_nodes, order_nodes.end(), [](const Node* a, const Node* b) {
        return a->order < b->order;
        });

    for (int i = offset_visible_order_nodes; i < order_nodes.size(); i++) {
        order_nodes[i]->order = i;
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
            for (auto& u_neighbour : order_nodes[u]->edges) {
                for (auto& v_neighbour : order_nodes[v]->edges) {
                    if (u_neighbour.neighbour_id > v_neighbour.neighbour_id) {
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
            for (auto& i : order_nodes[u]->edges) {
                for (auto& j : order_nodes[v]->edges) {
                    if (i.neighbour_id > j.neighbour_id) {
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
            for (int i = order_nodes[u]->offset_visible_nodes; i < order_nodes[u]->edges.size(); i++) {
                for (int j = order_nodes[v]->offset_visible_nodes; j < order_nodes[v]->edges.size(); j++) {
                    if (order_nodes[u]->edges[i].neighbour_id > order_nodes[v]->edges[j].neighbour_id) {
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
    for (auto& i : order_nodes[a]->edges) {
        for (auto& j : order_nodes[b]->edges) {
            if (i.neighbour_id > j.neighbour_id) {
                crossings++;
            }
        }
    }
    return crossings;
}

//beachte: berechnet crossings ohne die crossings in den twins selber
long Graph::countCrossingsWithEdgeWeights() {
    if (m == 0) {
        return 0;
    }
    long crossings = 0;
    //iterate through visible movable_nodes
    for (int u = offset_visible_order_nodes; u < order_nodes.size() - 1; u++) {
        //iterate through visible movable_nodes+1
        for (int v = u + 1; v < order_nodes.size(); v++) {
            for (auto& u_neighbour : order_nodes[u]->edges) {
                for (auto& v_neighbour : order_nodes[v]->edges) {
                    if (u_neighbour.neighbour_id > v_neighbour.neighbour_id) {
                        crossings += u_neighbour.edge_weight * v_neighbour.edge_weight;
                    }
                }
            }
        }
    }
    return crossings;
}

void Graph::sortNeighbours() {
    for (int i = 0; i < graph.size(); i++) {
        std::sort(graph[i].edges.begin(), graph[i].edges.end(), compareNeighbours);
    }
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
    assert(order_nodes[order_of_node]->offset_visible_nodes != order_nodes[order_of_node]->edges.size());

    Node* node = order_nodes[order_of_node];
    this->activeEdges -= node->edges.size();
    for (int i = node->offset_visible_nodes; i < node->edges.size(); i++) {
        Node* neighbour = &graph[node->edges[i].neighbour_id];

        //make node invisible in all adjacencie lists of fixed neighbours
        for (int j = neighbour->offset_visible_nodes; j < neighbour->edges.size(); j++) {
            if (neighbour->edges[j].neighbour_id != node->id) {
                continue;
            }
            //add node to invisible nodes and add visible counter
            std::swap(neighbour->edges[j].neighbour_id, neighbour->edges[neighbour->offset_visible_nodes].neighbour_id);
            neighbour->offset_visible_nodes++;
        }
    }

    //set visible counter to end of adjacency array -> as if all edges removed
    node->offset_visible_nodes = node->edges.size();

    //move invisible node to beginning of order_nodes
    while (node->order != this->offset_visible_order_nodes) {
        //swap with node one in front
        this->swapNodes(node->order, node->order - 1);
    }
    this->offset_visible_order_nodes++;
}

void Graph::makeNodeInvisibleBranching(int order_of_node) {
    Node* node = order_nodes[order_of_node];
    for (int i = node->offset_visible_nodes; i < node->edges.size(); i++) {
        Node* neighbour = &graph[node->edges[i].neighbour_id];

        //make node invisible in all adjacencie lists of fixed neighbours
        for (int j = neighbour->offset_visible_nodes; j < neighbour->edges.size(); j++) {
            if (neighbour->edges[j].neighbour_id != node->id) {
                continue;
            }
            //add node to invisible nodes and asdd visible counter
            std::swap(neighbour->edges[j].neighbour_id, neighbour->edges[neighbour->offset_visible_nodes].neighbour_id);
            neighbour->offset_visible_nodes++;
        }
    }

    //set visible counter to end of adjacency array -> as if all edges removed
    node->offset_visible_nodes = node->edges.size();
    graph[node->id].offset_visible_nodes = node->edges.size();
}

void Graph::makeNodeVisibleMarlon() {
    assert(offset_visible_order_nodes >= 1);

    int order_of_node = this->offset_visible_order_nodes - 1;
    //Make sure node is really invisible
    assert(order_nodes[order_of_node]->offset_visible_nodes == order_nodes[order_of_node]->edges.size());

    Node* node = order_nodes[order_of_node];
    activeEdges += node->edges.size();
    node->offset_visible_nodes = 0;

    for (int i = 0; i < node->edges.size(); i++) {
        Node* neighbour = &graph[node->edges[i].neighbour_id];

        //make node visible in all adjacencie lists of fixed neighbours
        for (int j = 0; j < neighbour->offset_visible_nodes; j++) {
            if (neighbour->edges[j].neighbour_id != node->id) {
                continue;
            }
            //add node to visible nodes and reduce visible counter
            std::swap(neighbour->edges[j].neighbour_id, neighbour->edges[neighbour->offset_visible_nodes - 1].neighbour_id);
            neighbour->offset_visible_nodes--;
        }
    }

    this->offset_visible_order_nodes--;
}

void Graph::makeNodeVisible(int order_of_node) {
    Node* node = order_nodes[order_of_node];
    node->offset_visible_nodes = 0;
    graph[node->id].offset_visible_nodes = 0;

    for (int i = 0; i < node->edges.size(); i++) {
        Node* neighbour = &graph[node->edges[i].neighbour_id];

        //make node visible in all adjacencie lists of fixed neighbours
        for (int j = 0; j < neighbour->offset_visible_nodes; j++) {
            if (neighbour->edges[j].neighbour_id != node->id) {
                continue;
            }
            //add node to visible nodes and reduce visible counter
            std::swap(neighbour->edges[j].neighbour_id, neighbour->edges[neighbour->offset_visible_nodes - 1].neighbour_id);
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

void Graph::MedianHeuristic() {
    // Iterate through moveable nodes
    for (int i = offset_visible_order_nodes; i < order_nodes.size(); i++) {
        auto& current_node = order_nodes[i];
        current_node->median_pos = 0;
        // If there are neighbours, select the median id
        if (current_node->offset_visible_nodes != current_node->edges.size()) {
            std::vector<int> neighbourIDs(0);
            for (int i = current_node->offset_visible_nodes; i < current_node->edges.size(); i++) {
                neighbourIDs.push_back(current_node->edges[i].neighbour_id);
            }

            int pos_index = ceil(neighbourIDs.size() / 2.0) - 1;
            current_node->median_pos = neighbourIDs[pos_index];
            // If even degree, further to the right
            if (neighbourIDs.size() % 2 == 0) {
                current_node->median_pos += 0.1;
            }
        }
    }

    std::sort(order_nodes.begin() + offset_visible_order_nodes, order_nodes.end(), [](const Node* a, const Node* b) {
        return a->median_pos < b->median_pos;
        });

    for (int i = offset_visible_order_nodes; i < order_nodes.size(); i++) {
        order_nodes[i]->order = i;
    }
}

void Graph::BarycenterHeuristicMarlon() {
    for (int i = offset_visible_order_nodes; i < order_nodes.size(); i++) {
        order_nodes[i]->barycenter_pos = 0;
        for (int j = order_nodes[i]->offset_visible_nodes; j < order_nodes[i]->edges.size(); j++) {
            order_nodes[i]->barycenter_pos += order_nodes[i]->edges[j].neighbour_id;
        }
        if ((order_nodes[i]->edges.size() - order_nodes[i]->offset_visible_nodes) != 0) {
            order_nodes[i]->barycenter_pos = order_nodes[i]->barycenter_pos / (order_nodes[i]->edges.size() - order_nodes[i]->offset_visible_nodes);
            //std::cout << "id: " << order_nodes[i]->id << " barycenter_pos: " << order_nodes[i]->barycenter_pos << std::endl;
        }
        else {
            order_nodes[i]->barycenter_pos = 0;
        }

    }

    std::sort(order_nodes.begin() + offset_visible_order_nodes, order_nodes.end(), [](const Node* a, const Node* b) {
        return a->barycenter_pos < b->barycenter_pos;
        });

    for (int i = offset_visible_order_nodes; i < order_nodes.size(); i++) {
        order_nodes[i]->order = i;
    }
}

void Graph::Barycenter_Heuristic() {
    for (int i = 0; i < order_nodes.size(); i++) {
        for (int j = 0; j < order_nodes[i]->edges.size(); j++) {
            order_nodes[i]->barycenter_pos += order_nodes[i]->edges[j].neighbour_id;
        }

        order_nodes[i]->barycenter_pos = order_nodes[i]->barycenter_pos / order_nodes[i]->edges.size();
        //std::cout << "id: " << order_nodes[i]->id << " barycenter_pos: " << order_nodes[i]->barycenter_pos << std::endl;
    }

    std::sort(order_nodes.begin(), order_nodes.end(), [](const Node* a, const Node* b) {
        return a->barycenter_pos < b->barycenter_pos;
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
        for (int j = 0; j < order_nodes[i]->edges.size(); j++) {
            int id_check = this->order_nodes[i]->id;
            if (this->order_nodes[i]->edges[j].neighbour_id != check.graph[id_check].edges[j].neighbour_id) { return false; }
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
void Graph::APUtil(int start_node_id, std::vector<bool>& visited, std::vector<int>& disc, std::vector<int>& low, int& time, int& parent, std::vector<bool>& isAP)
{
    // Count of children in DFS Tree
    int children = 0;

    // Mark the current node as visited
    visited[start_node_id] = true;

    // Initialize discovery time and low value
    disc[start_node_id] = low[start_node_id] = ++time;

    // Go through all vertices adjacent to this
    for (int i = graph[start_node_id].offset_visible_nodes; i < graph[start_node_id].edges.size(); i++) {
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (!visited[graph[start_node_id].edges[i].neighbour_id]) {
            children++;
            APUtil(graph[start_node_id].edges[i].neighbour_id, visited, disc, low, time, start_node_id, isAP);

            // Check if the subtree rooted with v has
            // a connection to one of the ancestors of u
            low[start_node_id] = std::min(low[start_node_id], low[graph[start_node_id].edges[i].neighbour_id]);

            // If start_node_id is not root and low value of one of
            // its child is more than discovery value of u.
            if (parent != -1 && low[graph[start_node_id].edges[i].neighbour_id] >= disc[start_node_id] && start_node_id < n0)
                isAP[start_node_id] = true;
        }

        // Update low value of start_node_id for parent function calls.
        else if (graph[start_node_id].edges[i].neighbour_id != parent)
            low[start_node_id] = std::min(low[start_node_id], disc[graph[start_node_id].edges[i].neighbour_id]);
    }

    // If start_node_id is root of DFS tree and has two or more children.
    if (parent == -1 && children > 1)
        isAP[start_node_id] = true;
}

bool compareTempPartition(const Partition_Intervall& a, const Partition_Intervall& b) {
    return a.interval_low < b.interval_low;
}


//change function such that it works for degree 0 nodes.
void Graph::DFS_AP_nodes(int start_node, std::vector<bool>& visited, std::vector<Partition_Intervall>& partition_intervall, std::vector<bool>& isAP) {
    // need to make sure at the beginning that all nodes with degree 0 are removed. we can do this since by setting their
    // order to be the right most position and adding an ignore marker.
    Stack stack;
    stack.push(start_node);

    std::vector<bool> AP_pushed_to_partition(n0, false);

    //get the index offset, to know where to where to start looping from to get the visible nodes.
    int index = graph[start_node].offset_visible_nodes;

    //Create the interval, set the high and low id to being the start_node.
    Partition_Intervall partition(graph[start_node].edges[index].neighbour_id, graph[start_node].edges[index].neighbour_id);

    //DFS
    while (!stack.isEmpty()) {
        int current_node = stack.peek();
        stack.pop();

        if (!visited[current_node]) {
            //if not visited add the visited node to the partition interval and set the visited to true.
            partition.nodes.push_back(&graph[current_node]);
            visited[current_node] = true;
        }

        //Explore neighbours of current node
        for (int i = graph[current_node].offset_visible_nodes; i < graph[current_node].edges.size(); i++) {
            int neighbour = graph[current_node].edges[i].neighbour_id;

            //update the interval high and low if necessary
            if (neighbour < n0 && neighbour < partition.interval_low) { partition.interval_low = neighbour; }
            if (neighbour < n0 && neighbour > partition.interval_high) { partition.interval_high = neighbour; }

            //if the neighbour is an AP only add it once to the partition vector 
            if (isAP[neighbour] && !AP_pushed_to_partition[neighbour]) {
                partition.nodes.push_back(&graph[neighbour]);
                AP_pushed_to_partition[neighbour] = true;
                continue;
            }

            if (!visited[neighbour]) {
                stack.push(neighbour);
            }
        }
    }
    std::sort(partition.nodes.begin(), partition.nodes.end(), compareNodeID);
    partition_intervall.push_back(partition);
}


void Graph::AP_Intervall() {
    //static int counter = 0;
    //std::cout << "new AP call " << counter++ << std::endl;

    std::vector<int> disc(this->n0 + this->n1, 0);
    std::vector<int> low(this->n0 + this->n1);
    std::vector<bool> visited(this->n0 + this->n1, false);
    std::vector<bool> isAP(this->n0 + this->n1, false);
    int time = 0, par = -1;

    partitions.clear();

    // Adding this loop so that the code works even if we are given disconnected graph
    for (int start_node_id = 0; start_node_id < this->n0; start_node_id++) {
        if (!visited[start_node_id]) {
            APUtil(start_node_id, visited, disc, low, time, par, isAP);
        }
    }

    visited = std::vector<bool>(this->n0 + this->n1, false); //need to use visited again for the interval
    std::vector<Partition_Intervall> temp_partitions;

    for (int node_ix = 0; node_ix < n0; node_ix++) {
        int node = graph[node_ix].id;
        for (int neighbour_ix = graph[node].offset_visible_nodes; neighbour_ix < graph[node].edges.size(); neighbour_ix++) {
            if (!visited[graph[node].edges[neighbour_ix].neighbour_id]) {
                //Run DFS add all nodes that are reachable from this neighbour node. Do not add AP marked nodes to stack. 
                // return the nodes reachable from the neighbour node as a vector, and return the interval of the fixed nodes. do so as a pair.
                // add the returned data to the temp_partitions
                DFS_AP_nodes(graph[node].edges[neighbour_ix].neighbour_id, visited, temp_partitions, isAP);
            }
        }
    }
    //sort temp partitions by last interval_high.
    std::sort(temp_partitions.begin(), temp_partitions.end(), compareTempPartition);

    //all neighbours have been visited + intervall has been found.
    for (int first_partition_ix = 0; first_partition_ix < temp_partitions.size(); first_partition_ix++) {
        //select the first partition to make the comparison with
        Partition_Intervall* first_partition = &temp_partitions[first_partition_ix];
        if (first_partition->ignore) { continue; }

        //select the second partition which is the next after the first partition
        for (int second_partition_ix = first_partition_ix + 1; second_partition_ix < temp_partitions.size(); second_partition_ix++) {
            Partition_Intervall* second_partition = &temp_partitions[second_partition_ix];
            //check if need to fuse
            if (!(first_partition->interval_high > second_partition->interval_low) || second_partition->ignore) { continue; }
            //fuse is needed, partitions are in the same interval
            //reset to 0 to start from the beginning.

            //CHANGE : could be that i need to add it.
            //first_partition_ix = 0;

            std::vector<Node*> first_partition_copy = first_partition->nodes;
            first_partition->nodes.clear();

            int first_index = 0;
            int second_index = 0;

            // Merge the two arrays
            while (first_index < first_partition_copy.size() && second_index < second_partition->nodes.size()) {
                if (first_partition_copy[first_index]->id == second_partition->nodes[second_index]->id) {
                    first_partition->nodes.push_back(first_partition_copy[first_index++]);
                    second_index++;
                }
                else if (first_partition_copy[first_index]->id < second_partition->nodes[second_index]->id) {
                    first_partition->nodes.push_back(first_partition_copy[first_index++]);
                }
                else {
                    first_partition->nodes.push_back(second_partition->nodes[second_index++]);
                }
            }

            while (first_index < first_partition_copy.size()) {
                first_partition->nodes.push_back(first_partition_copy[first_index++]);
            }
            while (second_index < second_partition->nodes.size()) {
                first_partition->nodes.push_back(second_partition->nodes[second_index++]);
            }
            
            //Change the interval's high and low if necessary after the merge
            if (first_partition->interval_low > second_partition->interval_low) { first_partition->interval_low = second_partition->interval_low; }
            if (first_partition->interval_high < second_partition->interval_high) { first_partition->interval_high = second_partition->interval_high; }
            second_partition->ignore = true;
        }
    }

    for (int partition_ix = 0; partition_ix < temp_partitions.size(); partition_ix++) {
        if (!temp_partitions[partition_ix].ignore) {
            partitions.push_back(temp_partitions[partition_ix].nodes);
        }
    }

}

void Graph::setNoPartitioning() {
    if (this->partitions.empty()) {
        std::vector<Node*> tmp(0);
        for (auto& node : this->graph) {
            tmp.push_back(&node);
        }
        this->partitions.push_back(tmp);
    }
}

void Graph::setOrderNodes(std::vector<Node*> order) {
    this->order_nodes = order;

    for (int i = offset_visible_order_nodes; i < order_nodes.size(); i++) {
        order_nodes[i]->order = i;
    }
}

void Branch_and_Bound(Graph* G) {
    G->Barycenter_Heuristic();
    Graph verifier = *G;
    std::vector<Node*>best_configuration = G->getOrderNodes();
    int best_solution = G->countCrossingsBranching();

    std::cout << "Median Heuristic : " << best_solution << std::endl;

    Graph G_branch = *G;
    for (int i = 1; i < G_branch.getOrderNodes().size(); i++) {
        G_branch.makeNodeInvisibleBranching(i);
    }

    exploreBranch(G_branch, *G, 0, best_solution, best_configuration);

    std::cout << "best solution " << best_solution << std::endl;
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

int calculateSpan(Node* node) {
    if (node->edges.empty()) {
        return 0;
    }

    // Find the minimum and maximum node IDs in the neighborhood
    auto minMaxIDs = std::minmax_element(node->edges.begin(), node->edges.end(),
        [](const Edge& a, const Edge& b) {
            return a.neighbour_id < b.neighbour_id;
        });

    int minID = minMaxIDs.first->neighbour_id;
    int maxID = minMaxIDs.second->neighbour_id;

    // Calculate the span
    return maxID - minID;
}

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
            for (int i = node->offset_visible_nodes; i < node->edges.size(); i++) {
                m++;
            }
        }
    }

    Graph* partGraph = new Graph(n0, n1, m); // use smartppointers

    for (int i = 0; i < partition.size(); i++) {
        partGraph->setOldID(i, partition[i]->id);
        //Need this so to not add forward AND backward edge
        if (partition[i]->id >= originalN0) {
            for (int j = partition[i]->offset_visible_nodes; j < partition[i]->edges.size(); j++) {
                partGraph->addEdge(i, std::distance(partitionIndices.begin(), std::find(partitionIndices.begin(), partitionIndices.end(), partition[i]->edges[j].neighbour_id)));
            }
        }
    }
    partGraph->sortNeighbours();
    return partGraph;
}

void TomAlvAlg(Graph& g) {
    
    // freenode position id -> pos
    // permutation pos -> id

    HeuristicGraph graphTomAlv = HeuristicGraph<int, int>(g.getN0(), g.getN1(), g.getM());

    for (int i = g.getOffsetVisibleOrderNodes(); i < g.getOrderNodes().size(); i++) {
        for (int j = 0; j < g.getOrderNodes()[i]->edges.size(); j++) {
            graphTomAlv.addEdge(g.getOrderNodes()[i]->id - g.getN0(), g.getOrderNodes()[i]->edges[j].neighbour_id);
        }
    }

    std::cout << "Crossings before TomAlv Algorithm : " << graphTomAlv.getCrossings() << std::endl;
    bool converged = heuristic_algorithm::HeuristicAlgorithm<HeuristicGraph<int, int>>(graphTomAlv, true, true, true);
    const std::vector<int>& permutation = graphTomAlv.getPermutation();
    std::cout << "Crossings After TomAlv Algorithm : " << graphTomAlv.getCrossings() << std::endl;

    
    std::vector<Node*> new_order(g.getOrderNodes().size());

    for (int i = 0; i < permutation.size(); i++) {
        new_order[i] = &g.getGraph()[permutation[i] + g.getN0()];
        new_order[i]->median_pos = i;
    }

    g.setOrderNodes(new_order);
    std::cout << "Crossing from our graph after TomAlv Algorithm : " << g.countCrossingsMarlon() << std::endl;

}

/*
std::pair<std::vector<Node*>, long> ExactSolution(Graph& g) {
    CrossGuard::Graph g_exact(g.getN0(), g.getN1());
    //std::cout << "Offset visible order nodes " << g.getOffsetVisibleOrderNodes() << std::endl;
    //std::cout << "Order Nodes Size : " << g.getOrderNodes().size() << std::endl;

    for (int i = g.getOffsetVisibleOrderNodes(); i < g.getOrderNodes().size(); i++) {
        for (int j = g.getOrderNodes()[i]->offset_visible_nodes; j < g.getOrderNodes()[i]->edges.size(); j++) {
            //TODO: Need to add Edge weights, done.
            
            //g_exact.add_edge(g.getOrderNodes()[i]->edges[j].neighbour_id, g.getOrderNodes()[i]->id, g.getOrderNodes()[i]->edges[j].edge_weight);
            //std::cout << "i : " << i << " j : " << j << std::endl;
            //std::cout << "i-te nodes amount of neighbours : " << g.getOrderNodes()[i]->edges.size() << std::endl;
            //std::cout << "i-te order node : " << g.getOrderNodes()[i]->id << std::endl;
            //std::cout << "j-te neighbour of i-te node: " << g.getOrderNodes()[i]->edges[j].neighbour_id << std::endl;

            g_exact.add_edge(g.getOrderNodes()[i]->edges[j].neighbour_id, g.getOrderNodes()[i]->id - g.getN0(), g.getOrderNodes()[i]->edges[j].edge_weight);
        }
    }

    //std::cout << "After adding edges" << std::endl;
    g_exact.finalize();

    CrossGuard::Solver s(g_exact);
    s.solve(true);
    CrossGuard::AlignedVector<CrossGuard::u32> solver_solution = s.get_solution();
    long sumCrossings = g_exact.determine_n_cuts(solver_solution);
    std::vector<Node*> new_order = g.getOrderNodes();
    
    for (int i = 0; i < solver_solution.size(); i++) {
        new_order[i] = &g.getGraph()[solver_solution[i] + g.getN0()];
    }

    std::cout << "Crossings with henning : " << sumCrossings << std::endl;
    std::cout << "Crossing from graph before new order : " << g.countCrossingsMarlon() << std::endl;
    g.setOrderNodes(new_order);
    std::cout << " Crossing from graph after new order: " << g.countCrossingsMarlon() << std::endl;

    if (sumCrossings != g.countCrossingsMarlon()) {
        std::cout << "SUM OF HENNING AND COUNTCROSSING FUNCTION DID NOT DELIVER THE SAME RESULT" << std::endl;
    }

    return std::make_pair(g.getOrderNodes(), g.countCrossingsMarlon());;
}
*/


std::pair<std::vector<Node*>, long> branching(Graph* g, std::vector<general_reduction*> reductionTypes, int method1, int method2, bool fast) {
    bool changed = false;
    int twins_count = 0;
    int almostTwins_count = 0;
    std::pair<std::vector<Node*>, long> result;

    for (auto& reduct : reductionTypes) {
        if (!g->getOptimal()) {
            if (reduct->get_reduction_type() == Twins) {
                twins_count = reduct->reduce(g);

                if (twins_count > 0){
                    changed = true;
                }
            } else if (reduct->get_reduction_type() == AlmostTwins){
                almostTwins_count = reduct->reduce(g);

                if (almostTwins_count > 0){
                    changed = true;
                }
            } else {
                changed = reduct->reduce(g);
            }
        }
    }

    if (g->getOptimal()) { //order already optimal, s.t. no recursion or bruteforce needed
        if (fast){
            result = std::make_pair(g->getOrderNodes(), 0);
        } else{
            result = std::make_pair(g->getOrderNodes(), g->countCrossingsMarlon());
        }
    }
    // TODO: Check this maybe with param
    else if (changed) {
        result = BranchAndReduce(g, reductionTypes, method1, method2, fast);
    }
    //Reduce our instance if no more reductions applicable
    else {
        // Randomize which method we use to remove node
        if (method1 == 3){
            // Seed the random number generator with the current time
            std::srand(std::time(0));

            // Generate a random number between 0 and 2
            method1 = std::rand() % 3;
        }

        auto& orderNodes = g->getOrderNodes();
        int visibleNodeOffset = g->getOffsetVisibleOrderNodes();

        if ((orderNodes.size() - visibleNodeOffset) > 2) {
            // TODO: Try param maybe?
            // Step 1: What kind of node do we want to remove?
            // Method 1A: Find highest degree visible moveable node
            Node* foundNode = nullptr;
            if (method1 == 0){
                foundNode = orderNodes[visibleNodeOffset];
                int maxDegree = orderNodes[visibleNodeOffset]->edges.size() - orderNodes[visibleNodeOffset]->offset_visible_nodes;
                for (int i = visibleNodeOffset + 1; i < orderNodes.size(); i++) {
                    if (maxDegree < orderNodes[i]->edges.size() - orderNodes[i]->offset_visible_nodes) {
                        maxDegree = orderNodes[i]->edges.size() - orderNodes[i]->offset_visible_nodes;
                        foundNode = orderNodes[i];
                    }
                }
            }
            // Method 1B: Find most spanning node
            else if (method1 == 1){
                foundNode = orderNodes[visibleNodeOffset];
                int maxSpan = calculateSpan(foundNode);
                for (int i = visibleNodeOffset + 1; i < orderNodes.size(); i++) {
                    int currentSpan = calculateSpan(orderNodes[i]);
                    if (maxSpan < currentSpan) {
                        maxSpan = currentSpan;
                        foundNode = orderNodes[i];
                    }
                }
            }
            //Method 1C: Find max combined value of span and degree
            else if (method1 == 2){
                foundNode = orderNodes[visibleNodeOffset];
                double maxCombined = (double)calculateSpan(foundNode) / (foundNode->edges.size() - foundNode->offset_visible_nodes);
                for (int i = visibleNodeOffset + 1; i < orderNodes.size(); i++) {
                    double currentCombined = (double)calculateSpan(orderNodes[i]) / (orderNodes[i]->edges.size() - orderNodes[i]->offset_visible_nodes);
                    if (maxCombined < currentCombined) {
                        maxCombined = currentCombined;
                        foundNode = orderNodes[i];
                    }
                }
            }
            
            // Remove node for now
            g->makeNodeInvisibleMarlon(foundNode->order);
            // Solve on remaining nodes
            result = BranchAndReduce(g, reductionTypes, method1, method2, fast);
            // Add node back
            g->makeNodeVisibleMarlon();

            //TODO: Try param here, so e.g. if less than 10 nodes left try every position (Method 2A)
            // Step 2: Where do we want to place the previously removed node?
            // Method 2A: Calculate Crossings for every possible position of node
            if (method2 == 0){
                int minCrossings = g->countCrossingsMarlon();
                auto minOrder = g->getOrderNodes();
                for (int i = visibleNodeOffset + 1; i < orderNodes.size(); i++) {
                    g->swapNodes(foundNode->order, i);
                    int crossings = g->countCrossingsMarlon();
                    if (crossings < minCrossings) {
                        minCrossings = crossings;
                        minOrder = g->getOrderNodes();
                    }
                }
                // Reconstruct min order
                g->setOrderNodes(minOrder);
                result = std::make_pair(g->getOrderNodes(), minCrossings);
            }
            // Method 2B: Place node at nearest median position
            else if (method2 == 1){
                orderNodes = g->getOrderNodes();        
                // Calculate median of node we want to place back
                int median = 0;
                if (foundNode->offset_visible_nodes != foundNode->edges.size()){
                    std::vector<int> neighbourIDs(0);
                    for (int i = foundNode->offset_visible_nodes; i < foundNode->edges.size(); i++){
                        neighbourIDs.push_back(foundNode->edges[i].neighbour_id);
                    }

                    int pos_index = ceil(neighbourIDs.size() / 2.0) - 1;
                    median = neighbourIDs[pos_index];
                    // If even degree, further to the right
                    if (neighbourIDs.size() % 2 == 0){
                        median += 0.1;
                    }
                }

                // Iterate through all other nodes to find closest median
                int minDifference = std::numeric_limits<int>::max();
                Node* closestNode = nullptr;
                for (int i = visibleNodeOffset + 1; i < orderNodes.size(); i++) {
                    auto& current_node = orderNodes[i];
                    current_node->median_pos = 0;
                    // If there are neighbours, select the median id
                    if (current_node->offset_visible_nodes != current_node->edges.size()){
                        std::vector<int> neighbourIDs(0);
                        for (int i = current_node->offset_visible_nodes; i < current_node->edges.size(); i++){
                            neighbourIDs.push_back(current_node->edges[i].neighbour_id);
                        }

                        int pos_index = ceil(neighbourIDs.size() / 2.0) - 1;
                        current_node->median_pos = neighbourIDs[pos_index];
                        // If even degree, further to the right
                        if (neighbourIDs.size() % 2 == 0){
                            current_node->median_pos += 0.1;
                        }

                        // Calculate the absolute difference between node median and the median of the node we want to place back
                        int difference = std::abs(current_node->median_pos - median);

                        // Update closest node if this node has a smaller difference
                        if (difference < minDifference) {
                            minDifference = difference;
                            closestNode = current_node;
                        }
                    }
                }
                foundNode->order = closestNode->order;
                g->sortOrderNodesByOrder();

                if (fast){
                    result = std::make_pair(g->getOrderNodes(), 0);
                } else {
                    result = std::make_pair(g->getOrderNodes(), g->countCrossingsMarlon());
                }
            // Method 2C: Calculate Crossings for every possible position of node using only crossings from that node
            } else if (method2 == 2){
                orderNodes = g->getOrderNodes();
                // Init for position = visibleNodeOffset
                int minCrossings = 0;
                auto minOrder = g->getOrderNodes();
                //iterate through visible movable_nodes != current node
                for (int v = visibleNodeOffset + 1; v < orderNodes.size(); v++) {
                    for (auto& current_neighbour : foundNode->edges) {
                        for (auto& v_neighbour : orderNodes[v]->edges) {
                            if (current_neighbour.neighbour_id > v_neighbour.neighbour_id) {
                                //minCrossings++;
                                minCrossings += current_neighbour.edge_weight * v_neighbour.edge_weight;
                            }
                        }
                    }
                }

                // Loop through every possible position
                for (int i = visibleNodeOffset + 1; i < orderNodes.size(); i++) {
                    g->swapNodes(foundNode->order, i);

                    int crossings = 0;
                    //iterate through visible movable_nodes != current node
                    for (int v = visibleNodeOffset; v < orderNodes.size(); v++) {
                        // v == current node, no crossings to itself
                        if (v == foundNode->order){
                            continue;
                        } 
                        // v to the left of current Node
                        else if (v < foundNode->order){
                            for (auto& current_neighbour : foundNode->edges) {
                                for (auto& v_neighbour : orderNodes[v]->edges) {
                                    if (current_neighbour.neighbour_id < v_neighbour.neighbour_id) {
                                        //crossings++;
                                        crossings += current_neighbour.edge_weight * v_neighbour.edge_weight;
                                    }
                                }
                            }
                        }
                        // v to the right of current Node
                        else if (v > foundNode->order){
                            for (auto& current_neighbour : foundNode->edges) {
                                for (auto& v_neighbour : orderNodes[v]->edges) {
                                    if (current_neighbour.neighbour_id > v_neighbour.neighbour_id) {
                                        //crossings++;
                                        crossings += current_neighbour.edge_weight * v_neighbour.edge_weight;
                                    }
                                }
                            }
                        }     
                    }

                    if (crossings < minCrossings) {
                        minCrossings = crossings;
                        minOrder = g->getOrderNodes();
                    }
                }
                // Reconstruct min order
                g->setOrderNodes(minOrder);
                if (fast){
                    result = std::make_pair(g->getOrderNodes(), 0);
                } else {
                    result = std::make_pair(g->getOrderNodes(), g->countCrossingsMarlon());
                }
            }
        }
        // For 2 remaining nodes, try which should come first
        else if ((orderNodes.size() - visibleNodeOffset) == 2) {
            int crossings1 = g->countCrossingsMarlon();
            g->swapNodes(visibleNodeOffset, visibleNodeOffset + 1);
            int crossings2 = g->countCrossingsMarlon();
            // go back to previous order if it was better
            if (crossings1 < crossings2) {
                g->swapNodes(visibleNodeOffset, visibleNodeOffset + 1);
                result = std::make_pair(g->getOrderNodes(), crossings1);
            }
            else {
                result = std::make_pair(g->getOrderNodes(), crossings2);
            }

        }
        // Less than two moveable nodes left -> nothing left to do
        else {
            result = std::make_pair(g->getOrderNodes(), 0);
        }
    }

    // Run through applies in backward order
    for (auto reduct = reductionTypes.rbegin(); reduct != reductionTypes.rend(); ++reduct) {
        bool applied = false;
        
        if ((*reduct)->get_reduction_type() == 3){
            applied = (*reduct)->apply(g, twins_count);
        } else if ((*reduct)->get_reduction_type() == 4){
            applied = (*reduct)->apply(g, almostTwins_count);
        }

        // Update result if there were changes made to order nodes
        if (applied) {
            result = std::make_pair(g->getOrderNodes(), g->countCrossingsMarlon());
        }
    }

    return result;
}

std::pair<std::vector<Node*>, long> BranchAndReduce(Graph* g, std::vector<general_reduction*> reductionTypes, int method1, int method2, bool fast) {
    //TODO: Try param here, maybe running Median once at beginning is often
    g->MedianHeuristic();
    //TomAlvAlg(*g);

    //Find partitions of Graph
    g->AP_Intervall();

    //g->printGraph();

    std::vector<std::vector<Node*>>& partitions = g->getPartitions();
    std::vector<Node*> solution(0);
    long sumCrossings = 0;
    // We could partition the graph
    if (partitions.size() > 1) {
        std::vector<std::pair<std::vector<Node*>, long>> results(0);

        //Get sub solutions
        for (auto& part : partitions) {
            Graph* partGraph = createGraphByPartition(g, part);

            //std::pair<std::vector<Node*>, long> result;

            //if (partGraph->getGraph().size() < 50) {
            //    result = ExactSolution(*partGraph);
            //}
            //else {
            //    result = branching(partGraph, reductionTypes, method1, method2, fast);
            //}

            auto result = branching(partGraph, reductionTypes, method1, method2, fast);
            results.push_back(result);
        }

        //Combine sub solutions
        std::vector<Node*> subSolutions(0);
        for (auto& result : results) {
            subSolutions.insert(subSolutions.end(), result.first.begin(), result.first.end());
            sumCrossings += result.second;
        }

        //Apply solution to original graph
        std::vector<Node*> oldOrder = g->getOrderNodes();
        auto& nodes = g->getGraph();

        //Fixed 0 edge issue: We remember all nodes that are invisible OR don't have any edges
        for (auto& node : oldOrder) {
            if (node->offset_visible_nodes == node->edges.size()) {
                solution.push_back(node);
            }
        }
        for (auto& node : subSolutions) {
            solution.push_back(&nodes[node->old_id]);
        }

        g->setOrderNodes(solution);

    }
    // We couldn't partition the graph
    else {

        //std::pair<std::vector<Node*>, long> result;

        //if (g->getGraph().size() < 50) {
        //    result = ExactSolution(*g);
        //}
        //else {
        //    result = branching(g, reductionTypes, method1, method2, fast);
        //}
        
        auto result = branching(g, reductionTypes, method1, method2, fast);
        solution = result.first;
        sumCrossings = result.second;

    }

    //TODO: Check if the solution is the order of the nodes.

    return std::make_pair(solution, sumCrossings);
}