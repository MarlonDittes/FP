#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

namespace median_algorithm {

// Here we don't copy the edges over, this spares some time
template <typename GraphType>
void medianAlgorithm(GraphType& myGraph) {
  using NT = typename GraphType::NodeType;

  std::vector<double> positions(myGraph.getFreeNodesSize(), 0.0);

  auto& edges = myGraph.getEdges();
  auto& permutationFreeNodes = myGraph.getFreeNodes();

  // Update positions and edges in-place
  for (NT i = 0; i < myGraph.getFreeNodesSize(); i++) {
    if (edges[i].size() == 0) {
      continue;
    }
    if (edges[i].size() % 2 == 0) {
      positions[i] =
          ((double)((edges[i][edges[i].size() / 2 - 1] + edges[i][edges[i].size() / 2 + 1]) / 2.0));
    } else {
      positions[i] = edges[i][edges[i].size() / 2];
    }
  }

  // Sort permutationFreeNodes based on positions
  std::sort(permutationFreeNodes.begin(), permutationFreeNodes.end(),
            [&positions](const NT& a, const NT& b) { return positions[a] < positions[b]; });
}
}  // namespace median_algorithm