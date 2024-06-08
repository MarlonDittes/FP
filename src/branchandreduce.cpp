#include "branchandreduce.h"
#include <memory>

/*
#include "../src_henning/src/definitions.h"
#include "../src_henning/src/macros.h"
#include "../src_henning/src/misc.h"
#include "../src_henning/src/graph_hen.h"
#include "../src_henning/src/solver_bf.h"
#include "../src_henning/src/exhaustive_solver.h"
#include "../src_henning/src/partitioner.h"
#include "../src_henning/src/solver.h"
#include "../src_henning/src/useless_reducer.h"
#include "../src_henning/src/front_back_reducer.h"
*/

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

std::unique_ptr<Graph> createGraphByPartition(Graph* g, std::vector<Node*> partition) {
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

    std::unique_ptr<Graph> partGraph = std::make_unique<Graph>(n0, n1, m);

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

    std::vector<int> mapping(g.getOrderNodes().size() - g.getOffsetVisibleOrderNodes());

    HeuristicGraph graphTomAlv = HeuristicGraph<int, int>(g.getN0(), g.getOrderNodes().size() - g.getOffsetVisibleOrderNodes(), g.getM());

    for (int i = g.getOffsetVisibleOrderNodes(); i < g.getOrderNodes().size(); i++) {
        g.getOrderNodes()[i]->TV_new_order_id = i - g.getOffsetVisibleOrderNodes();
        mapping[g.getOrderNodes()[i]->TV_new_order_id] = g.getOrderNodes()[i]->order;

        for (int j = g.getOrderNodes()[i]->offset_visible_nodes; j < g.getOrderNodes()[i]->edges.size(); j++) {
            graphTomAlv.addEdge(g.getOrderNodes()[i]->TV_new_order_id, g.getOrderNodes()[i]->edges[j].neighbour_id);
        }
    }

    bool converged = heuristic_algorithm::HeuristicAlgorithm<HeuristicGraph<int, int>>(graphTomAlv, true, true, true);
    const std::vector<int>& permutation = graphTomAlv.getPermutation();

    std::vector<Node*> new_order(g.getOrderNodes().size());

    for (int i = 0; i < g.getOffsetVisibleOrderNodes(); i++) {
        new_order[i] = g.getOrderNodes()[i];
        //new_order[i]->median_pos = i;
    }

    for (int i = 0; i < permutation.size(); i++) {
        new_order[i + g.getOffsetVisibleOrderNodes()] = g.getOrderNodes()[mapping[permutation[i]]];
        //new_order[i + g.getOffsetVisibleOrderNodes()]->median_pos = i;
    }

    g.setOrderNodes(new_order);

    for (int i = 0; i < g.getOrderNodes().size(); i++) {
        g.getOrderNodes()[i]->TV_new_order_id = -1;
    }

}

/*
std::pair<std::vector<Node*>, long> ExactSolution(Graph& g) {

    std::vector<int> mapping(g.getOrderNodes().size() - g.getOffsetVisibleOrderNodes());

    CrossGuard::Graph g_exact(g.getN0(), g.getOrderNodes().size() - g.getOffsetVisibleOrderNodes());

    for (int i = g.getOffsetVisibleOrderNodes(); i < g.getOrderNodes().size(); i++) {
        g.getOrderNodes()[i]->Hen_new_order_id = i - g.getOffsetVisibleOrderNodes();
        mapping[g.getOrderNodes()[i]->Hen_new_order_id] = g.getOrderNodes()[i]->order;

        for (int j = g.getOrderNodes()[i]->offset_visible_nodes; j < g.getOrderNodes()[i]->edges.size(); j++) {
            //TODO: Need to add Edge weights, done.

            g_exact.add_edge(g.getOrderNodes()[i]->edges[j].neighbour_id, g.getOrderNodes()[i]->Hen_new_order_id, g.getOrderNodes()[i]->edges[j].edge_weight);
        }
    }

    g_exact.finalize();

    CrossGuard::Solver s(g_exact);
    s.solve(false);
    CrossGuard::AlignedVector<CrossGuard::u32> solver_solution = s.get_solution();
    //long sumCrossings = g_exact.determine_n_cuts(solver_solution);

    std::vector<Node*> new_order = g.getOrderNodes();

    for (int i = 0; i < g.getOffsetVisibleOrderNodes(); i++) {
        new_order[i] = g.getOrderNodes()[i];
        //new_order[i]->median_pos = i;
    }

    for (int i = 0; i < solver_solution.size(); i++) {
        new_order[i + g.getOffsetVisibleOrderNodes()] = g.getOrderNodes()[mapping[solver_solution[i]]];
        //new_order[i + g.getOffsetVisibleOrderNodes()]->median_pos = i + g.getOffsetVisibleOrderNodes();
    }

    g.setOrderNodes(new_order);

    //if(sumCrossings != g.countCrossings()) {
    //    std::cout<<"CROSSINGS COUNT DO NOT MATCH"<<std::endl;
    //}

    for (int i = 0; i < g.getOrderNodes().size(); i++) {
        g.getOrderNodes()[i]->Hen_new_order_id = -1;
    }

    return std::make_pair(g.getOrderNodes(), g.countCrossings());

}
*/

int EXACT_SOLUTION_SIZE = 10;

std::pair<std::vector<Node*>, long> branching(Graph* g, std::vector<std::unique_ptr<general_reduction>>& reductionTypes, int method1, int method2, bool fast) {
    bool changed = false;
    int twins_count = 0;
    int almostTwins_count = 0;
    std::pair<std::vector<Node*>, long> result;

    for (auto& reduct : reductionTypes) {
        if (!g->getOptimal()) {
            if (reduct->get_reduction_type() == Twins) {
                twins_count = reduct->reduce(g);

                if (twins_count > 0) {
                    changed = true;
                }
            }
            else if (reduct->get_reduction_type() == AlmostTwins) {
                almostTwins_count = reduct->reduce(g);

                if (almostTwins_count > 0) {
                    changed = true;
                }
            }
            else {
                changed = reduct->reduce(g);
            }
        }
    }

    if (g->getOptimal()) { //order already optimal, s.t. no recursion or bruteforce needed
        if (fast) {
            result = std::make_pair(g->getOrderNodes(), 0);
        }
        else {
            result = std::make_pair(g->getOrderNodes(), g->countCrossings());
        }
    }
    // TODO: Check this maybe with param
    else if (changed) {
        result = BranchAndReduce(g, reductionTypes, method1, method2, fast);
    }
    //Reduce our instance if no more reductions applicable
    else {
        /*
        if (g->getOrderNodes().size() < EXACT_SOLUTION_SIZE) {
            result = ExactSolution(*g);
            return result;
        }*/



        // Randomize which method we use to remove node
        if (method1 == 3) {
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
            if (method1 == 0) {
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
            else if (method1 == 1) {
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
            else if (method1 == 2) {
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
            g->makeNodeInvisible(foundNode->order);
            // Solve on remaining nodes
            result = BranchAndReduce(g, reductionTypes, method1, method2, fast);
            // Add node back
            g->makeNodeVisible();

            //TODO: Try param here, so e.g. if less than 10 nodes left try every position (Method 2A)
            // Step 2: Where do we want to place the previously removed node?
            // Method 2A: Calculate Crossings for every possible position of node
            if (method2 == 0) {
                int minCrossings = g->countCrossings();
                auto minOrder = g->getOrderNodes();
                for (int i = visibleNodeOffset + 1; i < orderNodes.size(); i++) {
                    g->swapNodes(foundNode->order, i);
                    int crossings = g->countCrossings();
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
            else if (method2 == 1) {
                orderNodes = g->getOrderNodes();
                // Calculate median of node we want to place back
                int median = 0;
                if (foundNode->offset_visible_nodes != foundNode->edges.size()) {
                    std::vector<int> neighbourIDs(0);
                    for (int i = foundNode->offset_visible_nodes; i < foundNode->edges.size(); i++) {
                        neighbourIDs.push_back(foundNode->edges[i].neighbour_id);
                    }

                    int pos_index = ceil(neighbourIDs.size() / 2.0) - 1;
                    median = neighbourIDs[pos_index];
                    // If even degree, further to the right
                    if (neighbourIDs.size() % 2 == 0) {
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

                if (fast) {
                    result = std::make_pair(g->getOrderNodes(), 0);
                }
                else {
                    result = std::make_pair(g->getOrderNodes(), g->countCrossings());
                }
                // Method 2C: Calculate Crossings for every possible position of node using only crossings from that node
            }
            else if (method2 == 2) {
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
                        if (v == foundNode->order) {
                            continue;
                        }
                        // v to the left of current Node
                        else if (v < foundNode->order) {
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
                        else if (v > foundNode->order) {
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
                if (fast) {
                    result = std::make_pair(g->getOrderNodes(), 0);
                }
                else {
                    result = std::make_pair(g->getOrderNodes(), g->countCrossings());
                }
            }
        }
        // For 2 remaining nodes, try which should come first
        else if ((orderNodes.size() - visibleNodeOffset) == 2) {
            int crossings1 = g->countCrossings();
            g->swapNodes(visibleNodeOffset, visibleNodeOffset + 1);
            int crossings2 = g->countCrossings();
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

        if ((*reduct)->get_reduction_type() == 3) {
            applied = (*reduct)->apply(g, twins_count);
        }
        else if ((*reduct)->get_reduction_type() == 4) {
            applied = (*reduct)->apply(g, almostTwins_count);
        }

        // Update result if there were changes made to order nodes
        if (applied) {
            result = std::make_pair(g->getOrderNodes(), g->countCrossings());
        }
    }

    return result;
}

std::pair<std::vector<Node*>, long> BranchAndReduce(Graph* g, std::vector<std::unique_ptr<general_reduction>>& reductionTypes, int method1, int method2, bool fast) {
    //TODO: Try param here, maybe running Median once at beginning is often
    //g->MedianHeuristic();
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
            auto partGraph = createGraphByPartition(g, part);  // Take ownership of the new graph
            auto result = branching(partGraph.get(), reductionTypes, method1, method2, fast);
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

        auto result = branching(g, reductionTypes, method1, method2, fast);
        solution = result.first;
        sumCrossings = result.second;

    }

    return std::make_pair(solution, sumCrossings);
}