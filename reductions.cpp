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
    bool foundTwins = false;

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
                    foundTwins = true;
                    restore_vec.push_back({g->getNodeByOrder(i)->id, g->getNodeByOrder(j)->id}); //save twins
                    g->makeNodeInvisible(j); //j is made invisible 
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
    return foundTwins;
}