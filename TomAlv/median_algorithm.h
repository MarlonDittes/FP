#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
namespace median_algorithm {

    // Here we don't copy the edges over, this spares some time
    template <typename GraphType>
    void medianAlgorithm(GraphType& myGraph) {
        using NT = typename GraphType::NodeType;
        std::vector<int> positions(myGraph.getFreeNodesSize(), 0);
        auto& edges = myGraph.getEdges();
        auto& permutationFreeNodes = myGraph.getFreeNodes();
        // Update positions and edges in-place
        for (NT i = 0; i < myGraph.getFreeNodesSize(); i++) {
            if (edges[i].size() == 0) {
                continue;
            }
            int posTmp = ceil((edges[i][(edges[i].size() / 2)] - 1));
            positions[i] = posTmp;
        }

        // Sort permutationFreeNodes based on positions
        std::sort(permutationFreeNodes.begin(), permutationFreeNodes.end(),
            [&positions, &edges](const NT& a, const NT& b) {
                // If two free nodes have the same median, they are to be ordered according to their
                // average neighbour
                if (positions[a] == positions[b]) {
                    double averageA = 0;
                    double averageB = 0;
                    for (NT j = 0; j < edges[a].size(); ++j) {
                        averageA += edges[a][j];
                    }
                    averageA /= edges[a].size();
                    for (NT j = 0; j < edges[b].size(); ++j) {
                        averageB += edges[b][j];
                    }
                    averageB /= edges[b].size();

                    return averageA < averageB;
                }
                else {
                    return positions[a] < positions[b];
                }
            });
    }
}  // namespace median_algorithm