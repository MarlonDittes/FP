#include "reductions.h"

bool ZeroEdge_reduction::reduce(Graph* g) {
    if (g->getActiveEdges() == 0) {
        //std::cout << "Zero Edge works \n";
        g->setOptimalTrue();
        return true;
    }
    return false;
}

bool Complete_reduction::reduce(Graph* g) {
    int n0 = g->getN0();
    int n1 = g->getN1() - g->getOffsetVisibleOrderNodes();
    int m = g->getActiveEdges();

    if (n0 * n1 == m) {
        //std::cout << "Complete works \n";
        g->setOptimalTrue();
        return true;
    }
    return false;
}

bool ZeroCrossings_reduction::reduce(Graph* g) {
    int crossings = g->countCrossingsMarlon();
    if (crossings == 0) {
        //std::cout << "Zero Crossing works \n";
        g->setOptimalTrue();
        return true;
    }
    return false;
}

bool Twins_reduction::reduce(Graph* g) {
    bool found_Twins = false;

    for (int i = g->getOffsetVisibleOrderNodes(); i < g->getSizeOfOrder()-1; i++) {
        for (int j = i+1; j < g->getSizeOfOrder(); j++) {
            if (g->getNodeByOrder(i)->edges.size() == g->getNodeByOrder(i)->offset_visible_nodes ||
                g->getNodeByOrder(j)->edges.size() == g->getNodeByOrder(j)->offset_visible_nodes) {
                    continue;
            }
            //compare number of edges
            if (g->getNodeByOrder(i)->edges.size() != g->getNodeByOrder(j)->edges.size()) {
                continue;
            }
            //compare hash value
            if (g->getNodeByOrder(i)->hash != g->getNodeByOrder(j)->hash) {
                continue;
            }

            for (int k = 0; k < g->getNodeByOrder(i)->edges.size(); k++) {
                if (g->getNodeByOrder(i)->edges[k].neighbour_id != g->getNodeByOrder(j)->edges[k].neighbour_id) {
                    break;
                }

                //found twins
                if (k == g->getNodeByOrder(i)->edges.size()-1) {
                    //std::cout << "found Twins: " << g->getNodeByOrder(i)->id + 1 << " and "<< g->getNodeByOrder(j)->id + 1<< std::endl;
                    found_Twins = true;
                    restore_vec.push_back({g->getNodeByOrder(i)->id, g->getNodeByOrder(j)->id}); //save twins
                    for (auto& edge : g->getNodeByOrder(i)->edges) {
                        edge.edge_weight++;
                    }
                    g->makeNodeInvisibleMarlon(j);
                    std::cout << "g->getOffsetVisibleOrderNodes(): " << g->getOffsetVisibleOrderNodes() << std::endl;
                    i++;
                }
            }
        }
    }
    
    std::cout << "number of twins: " << restore_vec.size() << std::endl;
    for (int i = 0; i < restore_vec.size(); i++) {
        std::cout << "main: " << restore_vec[i].main << " , twin: " << restore_vec[i].twin << std::endl;
    }

    /*std::cout << "edge weights: " << std::endl;
    for (auto& node : g->getOrderNodes()) {
        std::cout << "node: " << node->id << ", ";
        for (auto& neighbour : node->neighbours) {
            std::cout << "neighbour: " << neighbour.neighbour_id << ", edge_weight: " << neighbour.edge_weight << std::endl;
        }
    }
    */
    return found_Twins;
}
void Twins_reduction::apply(Graph* g, int twinsCount){
    while (twinsCount > 0){
        g->makeNodeVisibleMarlon();
        auto pair = restore_vec[restore_vec.size()-1];
        restore_vec.pop_back();

        g->setOrderByNode(pair.twin, g->getOrderByNode(pair.main));
        twinsCount--;
    }
    g->sortOrderNodesByOrder();
}

bool AlmostTwin_reduction::reduce(Graph* g) {
    bool found_AlmostTwins = false;

    for (int i = g->getOffsetVisibleOrderNodes(); i < g->getSizeOfOrder()-1; i++) {
        for (int j = i+1; j < g->getSizeOfOrder(); j++) {
            if (g->getNodeByOrder(i)->offset_visible_nodes == g->getNodeByOrder(i)->edges.size() || 
                g->getNodeByOrder(j)->offset_visible_nodes == g->getNodeByOrder(j)->edges.size()) {
                    continue;
                }

            if (g->getNodeByOrder(i)->edges.size() == g->getNodeByOrder(j)->edges.size()+1
            || g->getNodeByOrder(i)->edges.size()+1 == g->getNodeByOrder(j)->edges.size()) {
                //find the node with less nodes
                int less = j; //node at order j has less neighbours -> i is main
                int more = i;
                if (g->getNodeByOrder(i)->edges.size() < g->getNodeByOrder(j)->edges.size()) {
                    less = i; //i has less nodes
                    more = j;
                }
                //first neighbour identical -> if j is almost twin it has to be on the left side of i
                if (g->getNodeByOrder(i)->edges[0].neighbour_id == g->getNodeByOrder(j)->edges[0].neighbour_id) {
                    // need tmp variable for this edge case
                    bool tmp = false;
                    if (g->getNodeByOrder(less)->edges.size() == 1) {
                        tmp = true;
                    }
                    else {
                        for (int k = 1; k < g->getNodeByOrder(less)->edges.size(); k++) {
                            if (g->getNodeByOrder(i)->edges[k].neighbour_id != g->getNodeByOrder(j)->edges[k].neighbour_id) { break; }
                            //found Almost Twins
                            if (k == g->getNodeByOrder(less)->edges.size()-1) { tmp = true; }
                        }
                    }

                    if (tmp) {
                        found_AlmostTwins = true;
                        restore_vec.push_back({g->getNodeByOrder(more)->id, g->getNodeByOrder(less)->id, 0});
                        //SET EDGE WEIGHTS  
                        for (int l = 0; l < g->getNodeByOrder(less)->edges.size(); l++) {
                            g->getNodeByOrder(more)->edges[l].edge_weight += g->getNodeByOrder(less)->edges[l].edge_weight;
                        }
                        g->makeNodeInvisibleMarlon(less);
                        i++;
                    }
                }
                //last neighbour identical -> if j almost twin, it has to be on the right of i 
                else if (g->getNodeByOrder(i)->edges[g->getNodeByOrder(i)->edges.size()-1].neighbour_id == g->getNodeByOrder(j)->edges[g->getNodeByOrder(j)->edges.size()-1].neighbour_id) {
                    bool tmp = false;
                    if (g->getNodeByOrder(less)->edges.size() == 1) {
                        tmp = true;
                    }
                    else {
                        //check if all other neighbours before are the same
                        for (int k = 0; k < g->getNodeByOrder(less)->edges.size()-1; k++) {
                            if (g->getNodeByOrder(less)->edges[k].neighbour_id != g->getNodeByOrder(more)->edges[k+1].neighbour_id) { break; }
                            //found Almost Twins
                            if (k == g->getNodeByOrder(less)->edges.size()-2) { tmp = true; }
                        }
                    }

                    if (tmp) {
                        found_AlmostTwins = true;
                        restore_vec.push_back({g->getNodeByOrder(more)->id, g->getNodeByOrder(less)->id, 1});
                        //SET EDGE WEIGHTS
                        for (int l = 1; l < g->getNodeByOrder(more)->edges.size(); l++) {
                            g->getNodeByOrder(more)->edges[l].edge_weight += g->getNodeByOrder(less)->edges[l-1].edge_weight;
                        }
                        g->makeNodeInvisibleMarlon(less);
                        i++;
                    }
                      
                }

            }
        }
    }

    /*std::cout << "number of almost_twins: " << restore_vec.size() << std::endl;
    for (int i = 0; i < restore_vec.size(); i++) {
        std::cout << "main: " << restore_vec[i].main << " , twin: " << restore_vec[i].twin << ", side: " << restore_vec[i].side << std::endl;
    }
    std::cout << "edge weights: " << std::endl;
    for (auto& node : g->getOrderNodes()) {
        std::cout << "node: " << node->id << ", ";
        for (auto& neighbour : node->edges) {
            std::cout << "neighbour: " << neighbour.neighbour_id << ", edge_weight: " << neighbour.edge_weight << std::endl;
        }
    }*/

    return found_AlmostTwins;
}

bool Domination_reduction::reduce(Graph* g) {
    bool found_Domination = false;
    int p = 2; //parameter

    while (p <= 10) {
        for (int i = g->getOffsetVisibleOrderNodes(); i < g->getSizeOfOrder()-1; i++) {
            for (int j = i+1; j < g->getSizeOfOrder(); j++) {
                if (g->getNodeByOrder(i)->offset_visible_nodes == g->getNodeByOrder(i)->edges.size() || 
                    g->getNodeByOrder(j)->offset_visible_nodes == g->getNodeByOrder(j)->edges.size()) {
                        continue;
                    }
                if (g->getNodeByOrder(i)->edges.size() == g->getNodeByOrder(j)->edges.size()+p ||
                    g->getNodeByOrder(i)->edges.size()+p == g->getNodeByOrder(j)->edges.size()) {
                    //find the node with less nodes
                    int less = j; //node at order j has less neighbours -> i is main
                    int more = i;
                    if (g->getNodeByOrder(i)->edges.size() < g->getNodeByOrder(j)->edges.size()) {
                        less = i; //i has less nodes
                        more = j;
                    }
                    //find first identical node
                    for (int k = 0; k <= g->getNodeByOrder(more)->edges.size() - g->getNodeByOrder(less)->edges.size(); k++) {
                        bool tmp = false;
                        if (g->getNodeByOrder(more)->edges[k].neighbour_id != g->getNodeByOrder(less)->edges[0].neighbour_id) {
                            continue;
                        }
                        if (g->getNodeByOrder(less)->edges.size() == 1) { tmp = true; }
                        //check if all other neighbours the same
                        for (int l = 1; l < g->getNodeByOrder(less)->edges.size(); l++) { 
                            if (g->getNodeByOrder(more)->edges[k+l].neighbour_id != g->getNodeByOrder(less)->edges[l].neighbour_id) {
                                break;
                            }
                           //found Domination
                            if (l == g->getNodeByOrder(less)->edges.size()-1) { tmp = true; }
                        }

                        if (tmp) {
                            found_Domination = true;
                            if (k <= g->getNodeByOrder(more)->edges.size()-k-g->getNodeByOrder(less)->edges.size()) {
                                //place invisible node on the left side
                                restore_vec.push_back({g->getNodeByOrder(more)->id, g->getNodeByOrder(less)->id, 0, (int)g->getNodeByOrder(less)->edges.size(), g->getNodeByOrder(more)->edges[k].neighbour_id});
                            } else {
                                //place it on the right side
                                restore_vec.push_back({g->getNodeByOrder(more)->id, g->getNodeByOrder(less)->id, 1, (int)g->getNodeByOrder(less)->edges.size(), g->getNodeByOrder(more)->edges[k].neighbour_id});
                                }
                                //SET EDGE WEIGHTS
                            for (int m = 0; m < g->getNodeByOrder(less)->edges.size(); m++) {
                                g->getNodeByOrder(more)->edges[k+m].edge_weight += g->getNodeByOrder(less)->edges[m].edge_weight;
                            }
                            g->makeNodeInvisibleMarlon(less);

                            //!!! i++ nur wenn j Schleife nicht am Ende, sonst Fehler.. evtl. auch bei twins;
                            if (j < g->getSizeOfOrder()-1) {
                                i++;
                            }
                            break;
                        }
                    }
                }
            }
        }
        //increase parameter
        p++;
    }

/* std::cout << "number of domination: " << restore_vec.size() << std::endl;
    for (int i = 0; i < restore_vec.size(); i++) {
        std::cout << "main: " << restore_vec[i].main << " , part: " << restore_vec[i].part << ", side: " << restore_vec[i].side << ", size: " << restore_vec[i].domination_size << ", start: " << restore_vec[i].start << std::endl;
    }
    std::cout << "edge weights: " << std::endl;
    for (auto& node : g->getOrderNodes()) {
        std::cout << "node: " << node->id << ", ";
        for (auto& neighbour : node->edges) {
            std::cout << "neighbour: " << neighbour.neighbour_id << ", edge_weight: " << neighbour.edge_weight << std::endl;
        }
    }
*/
    return found_Domination;   
}
