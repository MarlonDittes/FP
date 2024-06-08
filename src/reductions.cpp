#include "reductions.h"

int ZeroEdge_reduction::reduce(Graph* g) {
    if (g->getActiveEdges() == 0) {
        //std::cout << "Zero Edge works \n";
        g->setOptimalTrue();
        return 1;
    }
    return 0;
}

int Complete_reduction::reduce(Graph* g) {
    int n0 = g->getN0();
    int n1 = g->getN1() - g->getOffsetVisibleOrderNodes();
    int m = g->getActiveEdges();

    if (n0 * n1 == m) {
        //std::cout << "Complete works \n";
        g->setOptimalTrue();
        return 1;
    }
    return 0;
}

int ZeroCrossings_reduction::reduce(Graph* g) {
    int crossings = g->countCrossings();
    if (crossings == 0) {
        //std::cout << "Zero Crossing works \n";
        g->setOptimalTrue();
        return 1;
    }
    return 0;
}

int Twins_reduction::reduce(Graph* g) {
    int found_Twins = 0;
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
                    found_Twins++;
                    restore_vec.push_back({g->getNodeByOrder(i)->id, g->getNodeByOrder(j)->id}); //save twins
                    for (auto& edge : g->getNodeByOrder(i)->edges) {
                        edge.edge_weight++;
                    }
                    g->makeNodeInvisible(j);
                    i++;
                }
            }
        }
    }
    /*
    std::cout << "number of twins: " << restore_vec.size() << std::endl;
    for (int i = 0; i < restore_vec.size(); i++) {
        std::cout << "main: " << restore_vec[i].main << " , twin: " << restore_vec[i].twin << std::endl;
    }
    */

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

bool Twins_reduction::apply(Graph* g, int twins_count){
    if (twins_count > 0){
        while (twins_count > 0){
            g->makeNodeVisible();
            auto pair = restore_vec[restore_vec.size()-1];
            restore_vec.pop_back();

            g->setOrderByNode(pair.twin, g->getOrderByNode(pair.main));
            twins_count--;
            // reset edge weights
            Node* main = g->getNodeByOrder(g->getOrderByNode(pair.main));
            for (auto& edge : main->edges){
                edge.edge_weight--;
            }
            
        }
        g->sortOrderNodesByOrder();
        return true;
    } else {
        return false;
    }
}

int Domination_reduction::reduce(Graph* g) {
    int found_domination = 0;
    //copy and sort order array after number of edges
    std::vector<Node*> sortedNodes = g->getOrderNodes();
    int offset = g->getOffsetVisibleOrderNodes();
    std::sort(sortedNodes.begin() + offset, sortedNodes.end(), [](const Node* a, const Node* b) {
        return a->edges.size() < b->edges.size();
    }); 
    int offset_sortedNodes = 0; //offset to nodes that not yet in domination
    //compare biggest with smallest degree
    for (int i = g->getN1()-1; i >= g->getOffsetVisibleOrderNodes(); i--) {
        for (int j = offset_sortedNodes; j < i; j++) {
            Node* more = sortedNodes[i];
            Node* less = sortedNodes[j];
            //if more and less same degree, no domination
            if (more->edges.size() == less->edges.size()) { break; }
            //find first identical node
            for (int k = 0; k <= more->edges.size() - less->edges.size(); k++) {
                bool tmp = false;
                if (more->edges[k].neighbour_id != less->edges[0].neighbour_id) { continue; }
                if (less->edges.size() == 1) { tmp = true; }
                //check if all other neighbours the same
                for (int l = 1; l < less->edges.size(); l++) { 
                    //OR condition: check neighbourhood of less if all neighbours next to each other, s. t. no other edges can be inside (NOT TESTED)
                    if (more->edges[k+l].neighbour_id != less->edges[l].neighbour_id
                        || less->edges[l].neighbour_id != less->edges[l-1].neighbour_id+1) { break; }
                   //found Domination
                    if (l == less->edges.size()-1) { tmp = true; }                        
                }

                if (tmp) {
                found_domination++;
                if (k <= more->edges.size()-k-less->edges.size()) {
                    //place invisible node on the left side
                    restore_vec.push_back({more->id, less->id, 0, (int)less->edges.size(), less->edges[0].neighbour_id, std::make_pair(k,k+less->edges.size()-1)});
                } else {
                    //place it on the right side
                    restore_vec.push_back({more->id, less->id, 1, (int)less->edges.size(), less->edges[less->edges.size()-1].neighbour_id, std::make_pair(k,k+less->edges.size()-1)});
                    }
                    //SET EDGE WEIGHTS
                    for (int m = 0; m < less->edges.size(); m++) {
                        more->edges[k+m].edge_weight += less->edges[m].edge_weight;
                    }
                    g->makeNodeInvisible(less->order);
                    //make node invisible in sortedNodes as well:
                    std::swap(sortedNodes[j], sortedNodes[offset_sortedNodes]);
                    offset_sortedNodes++;
                    break;
                }
            }
        }
    }

    /*std::cout << "number of domination: " << restore_vec.size() << std::endl;
        for (int i = 0; i < restore_vec.size(); i++) {
            std::cout << "main: " << restore_vec[i].main << " , part: " << restore_vec[i].part << ", side: " << restore_vec[i].side << ", size: " << restore_vec[i].domination_size << ", start: " << restore_vec[i].start << ", intervall: " << restore_vec[i].intervall.first << "-" << restore_vec[i].intervall.second << std::endl;
        }
        std::cout << "edge weights: " << std::endl;
        for (auto& node : g->getOrderNodes()) {
            std::cout << "node: " << node->id << ", ";
            for (auto& neighbour : node->edges) {
                std::cout << "neighbour: " << neighbour.neighbour_id << ", edge_weight: " << neighbour.edge_weight << std::endl;
            }
        }*/

    return found_domination;
}

bool Domination_reduction::apply(Graph* g, int domination_count) {
    if (domination_count > 0) {
        double offset = 1.0 / (domination_count + 1);
        while (domination_count > 0){
            auto& tuple = restore_vec[restore_vec.size()-1];
            int main = tuple.main;
            std::vector<restore_data> dominations_left;
            std::vector<restore_data> dominations_right;
            //pop all dominations with same main
            while (tuple.main == main && domination_count > 0 && !restore_vec.empty()) {
                if (tuple.side == 0) {
                    dominations_left.push_back(tuple);
                } else {
                    dominations_right.push_back(tuple);
                }
                restore_vec.pop_back();
                g->makeNodeVisible();
                //reset Edge Weights
                Node* partNode = g->getNodeByOrder(g->getOrderByNode(tuple.part));
                Node* mainNode = g->getNodeByOrder(g->getOrderByNode(tuple.main));
                //start value different for left and right side
                int i = 0;
                for (int m = tuple.intervall.first; m <= tuple.intervall.second; m++) {
                    mainNode->edges[m].edge_weight -= partNode->edges[i].edge_weight;
                    i++;
                }
                domination_count--;
                tuple = restore_vec[restore_vec.size()-1];
            }
            if (!dominations_left.empty()) {
                std::sort(dominations_left.begin(), dominations_left.end(), [](const restore_data& a, restore_data& b) {
                    //if same start value, decide by domination size
                    if (a.start != b.start) {
                        return a.start < b.start;
                    } else {
                        return a.domination_size < b.domination_size;
                    }
                    
                });
                for (int i = 1; i <= dominations_left.size(); i++) {
                    tuple = dominations_left[dominations_left.size()-i];
                    double orderOfPart = g->getOrderByNode(tuple.main) - i*offset;
                    g->setOrderByNode(tuple.part, orderOfPart);
                }
            }
            if (!dominations_right.empty()) {
                std::sort(dominations_right.begin(), dominations_right.end(), [](const restore_data& a, restore_data& b) {
                    //if same start value, decide by domination size
                    if (a.start != b.start) {
                        return a.start < b.start;
                    } else {
                        //symmetric sorting 
                        return a.domination_size > b.domination_size;
                    }
                });
                for (int i = 1; i <= dominations_right.size(); i++) {
                    tuple = dominations_right[i-1];
                    double orderOfPart = g->getOrderByNode(tuple.main) + i*offset;
                    g->setOrderByNode(tuple.part, orderOfPart);
                }
            }
        }
        g->sortOrderNodesByOrder();
        return true;
    }
    return false;
}

/*int AlmostTwin_reduction::reduce(Graph* g) {
    int found_AlmostTwins = 0;

    std::vector<Node*> copy = g->getOrderNodes();
    int offset = g->getOffsetVisibleOrderNodes();
    //copy and sort order array after number of edges -> find right almost twins (deterministic)
    std::sort(copy.begin() + offset, copy.end(), [](const Node* a, const Node* b) {
        return a->edges.size() < b->edges.size();
    }); 

    for (int i = offset; i < copy.size()-1; i++) {
        for (int j = i+1; j < copy.size(); j++) {

            Node* less = copy[j]; 
            Node* more = copy[i];
            
            if (copy[i]->edges.size() < copy[j]->edges.size()) {
                less = copy[i];
                more = copy[j];
            }

            if (less->edges.size()+1 == more->edges.size()) {
                
                //first neighbour identical -> if j is almost twin it has to be on the left side of i
                if (less->edges[0].neighbour_id == more->edges[0].neighbour_id) {
                    // need tmp variable for this edge case
                    bool tmp = false;
                    if (less->edges.size() == 1) {
                        tmp = true;
                    }
                    else {
                        for (int k = 1; k < less->edges.size(); k++) {
                            if (less->edges[k].neighbour_id != more->edges[k].neighbour_id) { break; }
                            //found Almost Twins
                            if (k == less->edges.size()-1) { tmp = true; }
                        }
                    }

                    if (tmp) {
                        found_AlmostTwins++;
                        this->usage_count++;
                        restore_vec.push_back({more->id, less->id, 0});
                        //SET EDGE WEIGHTS  
                        for (int l = 0; l < less->edges.size(); l++) {
                            more->edges[l].edge_weight += less->edges[l].edge_weight;
                        }
                        g->makeNodeInvisible(less->order);
                        if (j < copy.size()-1) {
                            i++;
                        }
                    }
                }
                //last neighbour identical -> if j almost twin, it has to be on the right of i 
                else if (less->edges[less->edges.size()-1].neighbour_id == more->edges[more->edges.size()-1].neighbour_id) {
                    bool tmp = false;
                    if (less->edges.size() == 1) {
                        tmp = true;
                    }
                    else {
                        //check if all other neighbours before are the same
                        for (int k = 0; k < less->edges.size()-1; k++) {
                            if (less->edges[k].neighbour_id != more->edges[k+1].neighbour_id) { break; }
                            //found Almost Twins
                            if (k == less->edges.size()-2) { tmp = true; }
                        }
                    }

                    if (tmp) {
                        found_AlmostTwins++;
                        this->usage_count++;
                        restore_vec.push_back({more->id, less->id, 1});
                        //SET EDGE WEIGHTS
                        for (int l = 1; l < more->edges.size(); l++) {
                            more->edges[l].edge_weight += less->edges[l-1].edge_weight;
                        }
                        g->makeNodeInvisible(less->order);
                        if (j < copy.size()-1) {
                            i++;
                        }
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
        for (auto& neighbour : node->edges) {
            std::cout << "neighbour: " << neighbour.neighbour_id << ", edge_weight: " << neighbour.edge_weight << std::endl;
        }
    }

    return found_AlmostTwins;
}*/

/*bool AlmostTwin_reduction::apply(Graph* g, int twins_count) {
    if (twins_count > 0){
        double offset = 1.0 / (twins_count + 1);
        while (twins_count > 0){
            g->makeNodeVisible();
            auto tuple = restore_vec[restore_vec.size()-1];
            restore_vec.pop_back();

            // If on the left, offset to the left (negative)
            if (tuple.side == 0){
                offset = -offset;
            }

            double orderOfTwin = g->getOrderByNode(tuple.main) + offset;

            g->setOrderByNode(tuple.twin, g->getOrderByNode(tuple.main) + offset);
            twins_count--;
        }
        g->sortOrderNodesByOrder();
        return true;
    } else {
        return false;
    }
}*/