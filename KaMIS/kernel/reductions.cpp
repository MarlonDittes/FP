/******************************************************************************
* reductions.cpp
*****************************************************************************/


#include "reductions.h"
#include "branch_and_reduce_algorithm.h"
#include "data_structure/flow_graph.h"
#include "algorithms/push_relabel.h"

#include <utility>


typedef branch_and_reduce_algorithm::IS_status IS_status;

bool neighborhood_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	size_t oldn = status.remaining_nodes;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {
			NodeWeight neighbor_weights = 0;

			for (NodeID u : status.graph[v]) {
				neighbor_weights += status.weights[u];
			}

			if (status.weights[v] >= neighbor_weights) {
				br_alg->set(v, IS_status::included);
			}
		}
	}

	//std::cout << "neighbor redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

bool twin_reduction::reduce(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& neighbors = br_alg->buffers[0];
	auto& twin_candidates_set = br_alg->set_1;
	auto& tmp_set = br_alg->set_2;
	size_t oldn = status.remaining_nodes;

	NodeID twin;
	NodeWeight neighbors_weight;

	for (size_t v_idx = 0; v_idx < marker.current_size(); v_idx++) {
		NodeID v = marker.current_vertex(v_idx);

		if (status.node_status[v] == IS_status::not_set) {
			neighbors.clear();
			neighbors_weight = 0;

			for (NodeID neighbor : status.graph[v]) {
				neighbors.push_back(neighbor);
				neighbors_weight += status.weights[neighbor];
			}

			if (status.weights[v] >= neighbors_weight) {
				br_alg->set(v, IS_status::included);
				continue;
			}

			twin_candidates_set.clear();
			bool candidates_empty = true;

			for (NodeID neighbor : status.graph[neighbors[0]]) {
				if (neighbor != v && br_alg->deg(neighbor) == neighbors.size()) {
					twin_candidates_set.add(neighbor);
					candidates_empty = false;
					twin = neighbor;
				}
			}

			for (size_t i = 1; i < neighbors.size() && !candidates_empty; i++) {
				NodeID neighbor = neighbors[i];
				tmp_set.clear();
				candidates_empty = true;

				for (NodeID candidate : status.graph[neighbor]) {
					if (twin_candidates_set.get(candidate)) {
						tmp_set.add(candidate);
						candidates_empty = false;
						twin = candidate;
					}
				}

				std::swap(twin_candidates_set, tmp_set);
			}

			if (candidates_empty)
				continue;

			if (status.weights[v] + status.weights[twin] >= neighbors_weight) {
				br_alg->set(v, IS_status::included);
				br_alg->set(twin, IS_status::included);
			} else {
				fold(br_alg, v, twin);
			}
		}
	}

	//std::cout << "twin redu -> " << (oldn - status.remaining_nodes) << std::endl;

	return oldn != status.remaining_nodes;
}

void twin_reduction::fold(branch_and_reduce_algorithm* br_alg, NodeID main, NodeID twin) {
	auto& status = br_alg->status;

	restore_vec.push_back({ main, twin });

	br_alg->set(twin, IS_status::folded, true);
	status.weights[main] += status.weights[twin];

	status.folded_queue.push_back(get_reduction_type());

	br_alg->add_next_level_node(main);
	br_alg->add_next_level_neighborhood(main);
}

void twin_reduction::restore(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto& data = restore_vec.back();

	br_alg->unset(data.twin);
	status.weights[data.main] -= status.weights[data.twin];

	restore_vec.pop_back();
}

void twin_reduction::apply(branch_and_reduce_algorithm* br_alg) {
	auto& status = br_alg->status;
	auto main = restore_vec.back().main;
	auto twin = restore_vec.back().twin;

	restore(br_alg);

	if (status.node_status[main] == IS_status::included) {
		status.node_status[twin] = IS_status::included;
	} else {
		status.node_status[twin] = IS_status::excluded;
	}
}