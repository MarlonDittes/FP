#include "reductions.h"

bool ZeroEdge_reduction::reduce(Graph* g) {
    if (g->getActiveEdges() == 0) {
        std::cout << "Zero Edge works \n";
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
        std::cout << "Complete works \n";
        g->setOptimalTrue();
        return true;
    }
    return false;
}

bool ZeroCrossings_reduction::reduce(Graph* g) {
    int crossings = g->countCrossingsMarlon();
    if (crossings == 0) {
        std::cout << "Zero Crossing works \n";
        g->setOptimalTrue();
        return true;
    }
    return false;
}

bool Twins_reduction::reduce(Graph* g) {
    bool foundTwins = false;

    for (int i = 0; i < g->getN1()-1; i++) {
        for (int j = i+1; j < g->getN1(); j++) {
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
                    g->makeNodeInvisibleMarlon(j); //j is made invisible 
                    g->getNodeByOrder(i)->multiplier++;
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