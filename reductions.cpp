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
                if (g->getNodeByOrder(i)->neighbours[k].neighbour_id != g->getNodeByOrder(j)->neighbours[k].neighbour_id) {
                    break;
                }

                //found twins
                if (k == g->getNodeByOrder(i)->neighbours.size()-1) {
                    found_Twins = true;
                    restore_vec.push_back({g->getNodeByOrder(i)->id, g->getNodeByOrder(j)->id}); //save twins
                    for (auto& edge : g->getNodeByOrder(i)->neighbours) {
                        edge.edge_weight++;
                    }
                    g->makeNodeInvisibleMarlon(j); //j is made invisible
                }
            }
        }
    }
    std::cout << "number of twins: " << restore_vec.size() << std::endl;
    for (int i = 0; i < restore_vec.size(); i++) {
        std::cout << "main: " << restore_vec[i].main << " , twin: " << restore_vec[i].twin << std::endl;
    }
    std::cout << "edge weights: " << std::endl;
    for (auto& node : g->getOrderNodes()) {
        std::cout << "node: " << node->id << ", ";
        for (auto& neighbour : node->neighbours) {
            std::cout << "neighbour: " << neighbour.neighbour_id << ", edge_weight: " << neighbour.edge_weight << std::endl;
        }
    }
    return found_Twins;
}

//NOT TESTED YET
bool AlmostTwin_reduction::reduce(Graph* g) {
    bool found_AlmostTwins = false;

    for (int i = g->getOffsetVisibleOrderNodes(); i < g->getSizeOfOrder()-1; i++) {
        for (int j = i+1; j < g->getSizeOfOrder(); j++) {
            if (g->getNodeByOrder(i)->neighbours.size() == g->getNodeByOrder(j)->neighbours.size()+1
            || g->getNodeByOrder(i)->neighbours.size()+1 == g->getNodeByOrder(j)->neighbours.size()) {
                std::cout << "IF STATEMENT TRUE" << std::endl;
                //find the node with less nodes
                int less = j; //node at order j has less neighbours -> i is main
                int more = i;
                if (g->getNodeByOrder(i)->neighbours.size() < g->getNodeByOrder(j)->neighbours.size()) {
                    less = i; //i has less nodes
                    more = j;
                }
                std::cout << "LESS: " << g->getNodeByOrder(less)->id << std::endl;
                std::cout << "MORE: " << g->getNodeByOrder(more)->id << std::endl;
                //first neighbour identical -> if j is almost twin it has to be on the left side of i
                if (g->getNodeByOrder(i)->neighbours[0].neighbour_id == g->getNodeByOrder(j)->neighbours[0].neighbour_id) {
                    std::cout << "FIRST NEIGHBOUR IDENTICAL" << std::endl;
                    // need tmp variable for this edge case
                    bool tmp = false;
                    if (g->getNodeByOrder(less)->neighbours.size() == 1) {
                        tmp = true;
                    }
                    else {
                        for (int k = 1; k < g->getNodeByOrder(less)->neighbours.size(); k++) {
                            if (g->getNodeByOrder(i)->neighbours[k].neighbour_id != g->getNodeByOrder(j)->neighbours[k].neighbour_id) { break; }
                            //found Almost Twins
                            if (k == g->getNodeByOrder(less)->neighbours.size()-1) { tmp = true; }
                        }
                    }

                    if (tmp) {
                        std::cout << " FOUND ALMOST TWINS" << std::endl;
                        found_AlmostTwins = true;
                        restore_vec.push_back({g->getNodeByOrder(more)->id, g->getNodeByOrder(less)->id, 0});
                        //SET EDGE WEIGHTS  
                        for (int l = 0; l < g->getNodeByOrder(less)->neighbours.size(); l++) {
                            g->getNodeByOrder(more)->neighbours[l].edge_weight += g->getNodeByOrder(less)->neighbours[l].edge_weight;
                        }
                        g->makeNodeInvisibleMarlon(less);
                        break;
                    }
                }
                //last neighbour identical -> if j almost twin, it has to be on the right of i 
                else if (g->getNodeByOrder(i)->neighbours[g->getNodeByOrder(i)->neighbours.size()-1].neighbour_id == g->getNodeByOrder(j)->neighbours[g->getNodeByOrder(j)->neighbours.size()-1].neighbour_id) {
                    std::cout << "LAST NEIGHBOUR IDENTICAL" << std::endl;
                    bool tmp = false;
                    if (g->getNodeByOrder(less)->neighbours.size() == 1) {
                        tmp = true;
                    }
                    else {
                        //check if all other neighbours before are the same
                        for (int k = 0; k < g->getNodeByOrder(less)->neighbours.size()-1; k++) {
                            if (g->getNodeByOrder(less)->neighbours[k].neighbour_id != g->getNodeByOrder(more)->neighbours[k+1].neighbour_id) { break; }
                            //found Almost Twins
                            if (k == g->getNodeByOrder(less)->neighbours.size()-2) { tmp = true; }
                        }
                    }

                    if (tmp) {
                        std::cout << " FOUND ALMOST TWINS" << std::endl;
                        found_AlmostTwins = true;
                        restore_vec.push_back({g->getNodeByOrder(more)->id, g->getNodeByOrder(less)->id, 1});
                        //SET EDGE WEIGHTS
                        for (int l = 1; l < g->getNodeByOrder(more)->neighbours.size(); l++) {
                            g->getNodeByOrder(more)->neighbours[l].edge_weight += g->getNodeByOrder(less)->neighbours[l-1].edge_weight;
                        }
                        g->makeNodeInvisibleMarlon(less);
                        break;
                    }
                      
                }

            }
        }
    }

    std::cout << "number of almost_twins: " << restore_vec.size() << std::endl;
    for (int i = 0; i < restore_vec.size(); i++) {
        std::cout << "main: " << restore_vec[i].main << " , twin: " << restore_vec[i].twin << ", side: " << restore_vec[i].side << std::endl;
    }
    std::cout << "edge weights: " << std::endl;
    for (auto& node : g->getOrderNodes()) {
        std::cout << "node: " << node->id << ", ";
        for (auto& neighbour : node->neighbours) {
            std::cout << "neighbour: " << neighbour.neighbour_id << ", edge_weight: " << neighbour.edge_weight << std::endl;
        }
    }

    return found_AlmostTwins;
}
