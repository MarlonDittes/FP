#pragma once

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

template <typename NT, typename CCT>
class HeuristicGraph {
  // for each free node holds its neighbours
  std::vector<std::vector<NT>> freeNodes;
  // for each fixed node holds its neighbours
  std::vector<std::vector<NT>> fixedNodes;
  // free nodes current position in the permutation
  std::vector<NT> freeNodesPosition;
  // for each free node u, the first place hold the crossing sum with free node positioned to the
  // left of u and the second hold the sum to its right
  std::vector<std::array<CCT, 2>> leftRightCrossingSum;
  // hold free node id in its corrent position in the pemutation
  std::vector<NT> permutation;

 public:
  NT currentNumNodes() { return getFreeNodesSize() + getFixedNodesSize(); }
  NT currentNumEdges() { return 0; }
  using EdgeType = NT;
  using NodeType = NT;
  using WeightType = NT;
  using CrossingCountType = CCT;

  HeuristicGraph(const std::vector<std::vector<NodeType>>& freeNodes,
                 const std::vector<std::vector<NodeType>>& fixedNodes)
      : freeNodes(freeNodes),
        fixedNodes(fixedNodes),
        freeNodesPosition(freeNodes.size()),
        leftRightCrossingSum(freeNodes.size()),
        permutation(freeNodes.size()) {
    std::iota(freeNodesPosition.begin(), freeNodesPosition.end(), 0);
    std::iota(permutation.begin(), permutation.end(), 0);
    computeCrossingSums();
  }

  HeuristicGraph(NodeType numFixedNodes, NodeType numFreeNodes, CrossingCountType edgeNum)
      : freeNodes(numFreeNodes, std::vector<NodeType>(0)),
        fixedNodes(numFixedNodes, std::vector<NodeType>(0)),
        freeNodesPosition(freeNodes.size()),
        leftRightCrossingSum(freeNodes.size()),
        permutation(freeNodes.size()) {
    std::iota(freeNodesPosition.begin(), freeNodesPosition.end(), 0);
    std::iota(permutation.begin(), permutation.end(), 0);
    computeCrossingSums();
  }
  void addEdge(NodeType source,
               NodeType target) {  // where source is the freeNode and target is the fixedNode
    freeNodes[source].push_back(target);
    fixedNodes[target].push_back(source);
    return;
  }

  NodeType getFixedNodesSize() const { return fixedNodes.size(); }

  NodeType getFixedNodeNeighboursSize(NodeType nodeID) const { return fixedNodes[nodeID].size(); }

  const auto& getFixedNodeNeighbours(NodeType fixedNodeID) const { return fixedNodes[fixedNodeID]; }

  NodeType getFreeNodesSize() const { return freeNodes.size(); }

  NodeType getFreeNodeNeighboursSize(NodeType nodeID) const { return freeNodes[nodeID].size(); }

  const auto& getFreeNodeNeighbours(NodeType freeNodeID) const { return freeNodes[freeNodeID]; }

  CrossingCountType getLeftCrossings(NodeType freeNodeID) const {
    return leftRightCrossingSum[freeNodeID][0];
  }

  CrossingCountType getRightCrossings(NodeType freeNodeID) const {
    return leftRightCrossingSum[freeNodeID][1];
  }

  const auto& getPermutation() const { return permutation; }

  std::vector<NodeType>& getFreeNodes() { return permutation; }

  const auto& getFreeNodesPosition() const { return freeNodesPosition; }

  const auto& getEdges() const { return freeNodes; }

  /**
  This function returns the free node that is in index i of the perutation
  @param i The index of the permutation
*/
  const auto getPermutatuinAtIndex(NodeType i) const { return permutation[i]; }

  void setFreeNodes(const std::vector<NodeType>& newPermutation) {
    for (NodeType i = 0; i < newPermutation.size(); ++i) {
      freeNodesPosition[newPermutation[i]] = i;
      permutation[i] = newPermutation[i];
    }
  }

  // copmute the number of crossings created by the edges from two free nodes (u, v)
  // when u is to the left of v
  CrossingCountType const computeUVcrossing(NodeType u, NodeType v) {
    NodeType crossingSum = 0;
    for (const auto uNeighbour : freeNodes[u]) {
      for (const auto vNeighbour : freeNodes[v]) {
        crossingSum += vNeighbour < uNeighbour;
      }
    }
    return crossingSum;
  }

  CrossingCountType getCrossings() {
    CrossingCountType crossingCount = 0;
    for (size_t u = 0; u < freeNodes.size(); ++u) {
      for (size_t v = u + 1; v < freeNodes.size(); ++v) {
        crossingCount += computeUVcrossing(permutation[u], permutation[v]);
      }
    }
    return crossingCount;
  }

  /**
  This function switches the positions of two neighboring free nodes. Assumes u < v.
  @param u The first free node
  @param v The second free node
  @param isConditional switch only if reduce crossings
  @return True if the switch
*/
  bool switchNeighbours(NodeType u, NodeType v, bool isConditional) {
    // the number of crossings created by the edges from (u, v)
    NodeType uvSum = computeUVcrossing(u, v);
    // the number of crossings created by the edges from (v, u)
    NodeType vuSum = computeUVcrossing(v, u);
    if (isConditional) {
      if (uvSum <= vuSum) {
        return false;
      }
    }
    leftRightCrossingSum[u][1] -= uvSum;
    leftRightCrossingSum[v][1] += vuSum;
    leftRightCrossingSum[u][0] += vuSum;
    leftRightCrossingSum[v][0] -= uvSum;
    std::swap(permutation[freeNodesPosition[u]], permutation[freeNodesPosition[v]]);
    std::swap(freeNodesPosition[u], freeNodesPosition[v]);
    return true;
  }

  // copmute for each free node u the number of crossings created by the edges from nodes to the
  // left of u  and to its right
  void computeCrossingSums() {
    for (NodeType i = 0; i < freeNodes.size(); ++i) {
      for (NodeType j = i + 1; j < freeNodes.size(); ++j) {
        if (freeNodesPosition[i] > freeNodesPosition[j]) {
          CrossingCountType crossing = computeUVcrossing(j, i);
          leftRightCrossingSum[i][0] += crossing;
          leftRightCrossingSum[j][1] += crossing;
        } else {
          CrossingCountType crossing = computeUVcrossing(i, j);
          leftRightCrossingSum[i][1] += crossing;
          leftRightCrossingSum[j][0] += crossing;
        }
      }
    }
  }
};