#include <iostream>
#include "graph.h"
#include "io.h"
#include "stack.h"
#include "performance.h"
#include "reductions.h"
#include "ClpSimplex.hpp"

void LP(Graph* g) {

	int offset = g->getN0();

	ClpSimplex model;

	std::vector<std::vector<double>> M(g->getN1(), std::vector<double>(g->getN1(), 0.0));
	for (int i = 0; i < g->getOrderNodes().size(); i++) {
		for (int j = i + 1; j < g->getOrderNodes().size(); j++) {
			M[g->getNodeByOrder(i)->id - offset][g->getNodeByOrder(j)->id - offset] = 1;
		}
	}

	std::vector<std::vector<int>> crosses(g->getN1(), std::vector<int>(g->getN1(), 0));
	for (int i = 0; i < g->getN1(); i++) {
		for (int j = 0; j < g->getN1(); j++) {
			if (i == j) {
				continue;
			}
			else if (i < j) {
				crosses[i][j] = g->countCrossingsMarlon();
			}
			else {
				g->swapNodes(g->getNodeByOrder(i)->id - offset, g->getNodeByOrder(j)->id - offset);
				crosses[i][j] = g->countCrossingsMarlon();
				g->swapNodes(g->getNodeByOrder(i)->id - offset, g->getNodeByOrder(j)->id - offset);
			}
		}
	}

	int totalVariables = g->getN1() * (g->getN1() - 1);
	model.resize(0, totalVariables);


	for (int i = 0; i < g->getN1(); ++i) {
		for (int j = 0; j < g->getN1(); ++j) {
			if (i != j) {
				int index = i * (g->getN1() - 1) + j - (j > i ? 1 : 0);
				model.setObjectiveCoefficient(index, crosses[i][j]);
				model.setColumnBounds(index, 0.0, 1.0);
			}
		}
	}

	for (int i = 0; i < g->getN1(); ++i) {
		for (int j = 0; j < g->getN1(); ++j) {
			if (i != j) {
				int index_ij = i * (g->getN1() - 1) + j - (j > i ? 1 : 0);
				for (int k = 0; k < g->getN1(); ++k) {
					if (k != i && k != j) {
						int index_jk = j * (g->getN1() - 1) + k - (k > j ? 1 : 0);
						int index_ik = i * (g->getN1() - 1) + k - (k > i ? 1 : 0);

						// Set up the row for the constraint m[i,j] + m[j,k] - m[i,k] <= 1 --> Transitivity
						int indices[] = { index_ij, index_jk, index_ik };
						double coefficients[] = { 1.0, 1.0, -1.0 };
						model.addRow(3, indices, coefficients, -COIN_DBL_MAX, 1.0);
					}
				}
			}
		}
	}

	for (int i = 0; i < g->getN1(); ++i) {
		for (int j = 0; j < g->getN1(); ++j) {
			if (i != j) {
				int index_ij = i * (g->getN1() - 1) + j - (j > i ? 1 : 0);
				int index_ji = j * (g->getN1() - 1) + i - (i > j ? 1 : 0);
				// m[i,j] + m[j,i] = 1
				int indices[] = { index_ij, index_ji };
				double coefficients[] = { 1.0, 1.0 };
				model.addRow(2, indices, coefficients, 1.0, 1.0);

			}
		}
	}

	model.primal();
	const double* solution = model.primalColumnSolution();

	std::vector<int> new_node_order(g->getN1());

	//for (int i = 0; i < g->getN1(); ++i) {
	//	for (int j = 0; j < g->getN1(); ++j) {
	//		int index = i * (g->getN1() - 1) + j - (j > i ? 1 : 0);
	//		std::cout << "m[" << i << "][" << j << "] = " << solution[index] << " ";
	//	}
	//	std::cout << std::endl;
	//} 

	// TODO: Indexing could still be a problem have a look at that. Sometimes it gives me for i == j 0 and sometimes 1,
	// Assume it has something to do with the indexing.

	//Print order of Moveable nodes

	for (int i = 0; i < g->getN1(); i++) {
		int sum_order = 0;
		for (int j = 0; j < g->getN1(); j++) {
			if (i == j) { continue; }
			int index = j * (g->getN1() - 1) + i - (i > j ? 1 : 0);
			if (solution[index] >= 0.5) {
				sum_order += 1;
			}
		}
		//std::cout << sum_order << std::endl;
		new_node_order[g->getNodeByOrder(i)->id - g->getN0()] = sum_order;
	}

	//for (int i = 0; i < g->getOrderNodes().size(); i++) {
	//	std::cout << new_node_order[i] << std::endl;
	//}

	std::vector<Node*> new_node_array_order(g->getN1());

	for (int i = 0; i < g->getN1(); i++) {
		new_node_array_order[new_node_order[i]] = g->getNodeByOrder(i);
	}

	//for (int i = 0; i < new_node_array_order.size(); i++) {
	//	std::cout << new_node_array_order[i]->id << std::endl;
	//}

	g->setOrderNodes(new_node_array_order);


	std::cout << "Crossings with new LP order : " << g->countCrossings() << std::endl;

}