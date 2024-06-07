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
    for (int i = 0; i < this->partitions.size(); i++) {
        std::cout << "partition : " << i << " node in partition : ";
        for (int j = 0; j < this->partitions[i].size(); j++) {
            std::cout << this->partitions[i][j]->id + 1 << " ";
        }
        std::cout << std::endl;
    }
}

/*void Graph::printGraphByPartitions() {
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
}*/

void Graph::setOrderNodes(std::vector<Node*> order) {
    this->order_nodes = order;

    for (int i = offset_visible_order_nodes; i < order_nodes.size(); i++) {
        order_nodes[i]->order = i;
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

void Graph::sortNeighbours() {
    for (int i = 0; i < graph.size(); i++) {
        std::sort(graph[i].edges.begin(), graph[i].edges.end(), compareNeighbours);
    }
}

long Graph::countCrossings() {
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

/*long Graph::countCrossingsWithEdgeWeights() {
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
}*/

void Graph::swapNodes(int node0, int node1) {
    std::swap(order_nodes[node0]->order, order_nodes[node1]->order);
    std::swap(order_nodes[node0], order_nodes[node1]);
}

void Graph::makeNodeInvisible(int order_of_node) {
    assert(offset_visible_order_nodes <= order_of_node);
    assert(order_of_node <= order_nodes.size() - 1);

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

void Graph::makeNodeVisible() {
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

/*bool compareNodeOrder(Node* a, Node* b) {
    return a->order < b->order;
}*/

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

/*std::pair<std::vector<Node*>, long> Graph::Greedy() {
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
}*/

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

/*std::pair<std::vector<Node*>, long> bruteForce(Graph* g) {
    std::vector<Node*> baseOrder = g->getOrderNodes();
    std::sort(baseOrder.begin(), baseOrder.end(), compareNodeID);

    std::vector<Node*> bestOrder = baseOrder;
    long bestCrossings = g->countCrossings();

    do {
        g->setOrderNodes(baseOrder);
        // Display the current permutation
        std::cout << "Permutation:";
        for (const Node* node : baseOrder) {
            std::cout << " order_nodes: " << node->order << ", ID: " << node->id;
        }
        std::cout << "\n";

        long crossings = g->countCrossings();
        //std::cout << "Crossings: " << crossings << std::endl;

        if (crossings < bestCrossings) {
            bestOrder = baseOrder;
            bestCrossings = crossings;
        }

    } while (std::next_permutation(baseOrder.begin(), baseOrder.end(), compareNodeID));

    return std::make_pair(bestOrder, bestCrossings);
}*/

/*std::pair<std::vector<Node*>, long> bruteForceOnSubgraph(Graph* g, int begin, int end) {
    std::vector<Node*> baseOrder = g->getOrderNodes();
    std::sort(baseOrder.begin() + begin, baseOrder.begin() + end, compareNodeID);

    std::vector<Node*> bestOrder = baseOrder;
    long bestCrossings = g->countCrossings();

    do {
        g->setOrderNodes(baseOrder);
        // Display the current permutation
        //std::cout << "Permutation:";
        for (const Node* node : baseOrder) {
            std::cout << " order_nodes: " << node->order << ", ID: " << node->id;
        }
        std::cout << "\n";

        long crossings = g->countCrossings();
        //std::cout << "Crossings: " << crossings << std::endl;

        if (crossings < bestCrossings) {
            bestOrder = baseOrder;
            bestCrossings = crossings;
        }

    } while (std::next_permutation(baseOrder.begin() + begin, baseOrder.begin() + end, compareNodeID));

    return std::make_pair(bestOrder, bestCrossings);
}*/