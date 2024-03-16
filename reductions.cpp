#include "reductions.h"

bool ZeroEdge_reduction::reduce(Graph* g) {
    if (g->getM() == 0) {
        g->setOptimalTrue();
        return true;
    }
    return false;
}

bool Complete_reduction::reduce(Graph* g) {
    int n0 = g->getN0();
    int n1 = g->getN1();
    int m = g->getM();

    if (n0 * n1 == m) {
        g->setOptimalTrue();
        std::cout << "GRAPH IS COMPLETE" << std::endl;
        return true;
    }
    return false;
}

bool ZeroCrossings_reduction::reduce(Graph* g) {
    int crossings = g->countCrossings();
    if (crossings == 0) {
        g->setOptimalTrue();
        return true;
    }
    return false;
}

bool Twins_reduction::reduce(Graph* g) {
    bool found_Twins = false;

    for (int i = 0; i < g->getSizeOfOrder()-1; i++) {
        for (int j = i+1; j < g->getSizeOfOrder(); j++) {
            if (g->getNodeByOrder(i)->neighbours.size() != g->getNodeByOrder(j)->neighbours.size()) {
                continue;
            }

            for (int k = g->getNodeByOrder(i)->offset_visible_nodes; k < g->getNodeByOrder(i)->neighbours.size(); k++) {
                if (g->getNodeByOrder(i)->neighbours[k] != g->getNodeByOrder(j)->neighbours[k]) {
                    break;
                }

                //found twins
                if (k == g->getNodeByOrder(i)->neighbours.size()-1) {
                    found_Twins = true;
                    restore_vec.push_back({g->getNodeByOrder(i)->id, g->getNodeByOrder(j)->id}); //save twins
                    g->makeNodeInvisible(j); //j is made invisible
                    //!!! CHANGE to set edge weights 
                    g->getNodeByOrder(i)->multiplier++; 
                    //!!! save crossings in i for both j and i, as twins i and j are reduced to one node
                }
            }
        }
    }
    std::cout << "number of twins: " << restore_vec.size() << std::endl;
    for (int i = 0; i < restore_vec.size(); i++) {
        std::cout << "main: " << restore_vec[i].main << " , twin: " << restore_vec[i].twin << std::endl;
    }
    return found_Twins;
}

//NOT TESTED YET
bool AlmostTwin_reduction::reduce(Graph* g) {
    bool found_AlmostTwins = false;
    //work as if offset order array is included
    int orderOffset;

    for (int i = orderOffset; i < g->getSizeOfOrder()-1; i++) {
        for (int j = i+1; j < g->getSizeOfOrder(); j++) {
            if (g->getNodeByOrder(i)->neighbours.size() == g->getNodeByOrder(j)->neighbours.size()+1
            || g->getNodeByOrder(i)->neighbours.size()+1 == g->getNodeByOrder(j)->neighbours.size()) {
                //find the node with less nodes
                int less = j; //node at order j has less neighbours -> i is main
                int more = i;
                if (g->getNodeByOrder(i)->neighbours.size() < g->getNodeByOrder(j)->neighbours.size()) {
                    less = i; //i has less nodes
                    more = j;
                }
                //first neighbour identical -> if j is almost twin it has to be on the left side of i
                if (g->getNodeByOrder(i)->neighbours[0] == g->getNodeByOrder(j)->neighbours[0]) {
                    //check if all other neighbours are the same
                    for (int k = 1; k < g->getNodeByOrder(less)->neighbours.size(); k++) {
                        if (g->getNodeByOrder(i)->neighbours[k] != g->getNodeByOrder(j)->neighbours[k]) {
                            break;
                        }
                        //found Almost Twins
                        if (k == g->getNodeByOrder(less)->neighbours.size()-1) {
                            found_AlmostTwins = true;
                            restore_vec.push_back({g->getNodeByOrder(more)->id, g->getNodeByOrder(less)->id, 0});
                            //SET EDGE WEIGHTS  
                        }
                    }
                }
                //last neighbour identical -> if j almost twin, it has to be on the right of i 
                else if (g->getNodeByOrder(i)->neighbours[g->getNodeByOrder(i)->neighbours.size()-1] == g->getNodeByOrder(j)->neighbours[g->getNodeByOrder(i)->neighbours.size()-1]) {
                    //check if all others neighbours before are the same
                    for (int k = 1; k < g->getNodeByOrder(less)->neighbours.size()-1; k++) {
                        if (g->getNodeByOrder(i)->neighbours[k] != g->getNodeByOrder(j)->neighbours[k]) {
                            break;
                        }
                        //found Almost Twins
                        if (k == g->getNodeByOrder(less)->neighbours.size()-2 || g->getNodeByOrder(less)->neighbours.size() == 1) {
                            found_AlmostTwins = true;
                            restore_vec.push_back({g->getNodeByOrder(more)->id, g->getNodeByOrder(less)->id, 1});
                            //SET EDGE WEIGHTS
                        }
                    }
                }
            }
        }
    }
    return found_AlmostTwins;
}

std::pair<std::vector<Node*>, long> BranchAndReduce(Graph* g, std::vector<general_reduction*> reductionTypes) {
    g->Median_Heuristic();

    g->Partition();
    g->AP();

    std::vector<std::vector<Node*>> partitions = g->getPartitions();           // NEED TO ADJUST THIS PROBABLY TO BE VECTOR<NODE*> NOT VECTOR<INT>!!, wait for shai AP
    std::vector<std::pair<std::vector<Node*>, long>> results(0);
    for (auto part : partitions) {
        Graph* partGraph = createGraphByPartition(g, part);

        bool changed = false;
        for (auto& reduct : reductionTypes) {
            changed = reduct->reduce(partGraph);
        }

        //Check if need to brute force since no more reductions were applicable
        std::pair<std::vector<Node*>, long> result;

        if (partGraph->getOptimal()) { //order already optimal, s.t. no recursion or bruteforce needed
            result = std::make_pair(part, partGraph->countCrossings()); //BETTER: save crossing in partGraph instead of calculating it again, since it was calculated in reductions already
        }
        else {
            //decide if we should do this:
            int numberOfVertices = g->getN0() + g->getN1();
            bool otherCondition = numberOfVertices <= BRUTE_CUTOFF;

            if (changed == false || otherCondition) {
                result = bruteForce(g);
            }
            else {
                result = BranchAndReduce(partGraph, reductionTypes);
            }  
        }
        
        results.push_back(result);
    }

    std::vector<Node*> solution(0);
    long sumCrossings = 0;
    for (auto result : results) {
        solution.insert(solution.end(), result.first.begin(), result.first.end());
        sumCrossings += result.second;
    }

    return std::make_pair(solution, sumCrossings);
}