#ifndef PACE2024EXACT_TWIN_REDUCER_H
#define PACE2024EXACT_TWIN_REDUCER_H

#include <numeric>
#include <ranges>
#include <utility>

#include "definitions.h"
#include "macros.h"
#include "misc.h"
#include "graph_hen.h"
#include "translation_table.h"

namespace CrossGuard {

    /**
     * Uses twin reductions to decrease the size of a graph.
     */
    class TwinReducer {
    private:
        const Graph &m_graph;

        // Variables for Twins
        AlignedVector<AlignedVector<u32>> m_twins;
        TranslationTable m_tt_twins;

    public:
        /**
         * Default constructor.
         *
         * @param g The graph.
         */
        explicit TwinReducer(const Graph &g) : m_graph(g) {}

        /**
         * Reduces the graph and returns a new graph.
         *
         * @return The newly reduced graph.
         */
        inline Graph reduce() {
            Graph g = m_graph;

            find_twins(g);
            g = reduce_twins(g);

            return g;
        }

        /**
         * The function will transform the solution of the reduced graph to the
         * unreduced graph.
         *
         * @param sol The solution of the reduced graph.
         * @return The solution to the original graph.
         */
        inline AlignedVector<u32> back_propagate(const AlignedVector<u32> &sol) {
            AlignedVector<u32> new_sol = sol;

            new_sol = back_propagate_twins(new_sol);

            return new_sol;
        }

    private:
        /**
         * Looks for sets of movable vertices that all have the same neighborhood.
         * Note that we name them twins, but in reality they can be of any size.
         * These twins can be placed one after the other and we only need to hold
         * one representative.
         */
        inline void find_twins(const Graph &g) {
            ASSERT(g.is_finalized);

            // collect hashes with vertices
            AlignedVector<std::pair<u64, u32>> hash_b(g.n_B);
            for(u32 i = 0; i < g.n_B; ++i){
                u32 vertex_b = i;
                u64 hash = g.adj_hash[i];
                hash_b[i] = {hash, vertex_b};
            }

            // sort the hashes
            std::sort(hash_b.begin(), hash_b.end(), [](auto &left, auto &right) {
                return left.first < right.first;
            });

            AlignedVector<u32> currentGroup;
            currentGroup.push_back(hash_b[0].second);

            for (size_t i = 1; i < g.n_B; ++i) {
                if (hash_b[i].first == hash_b[i - 1].first) {
                    // If the current element is equal to the previous one, add it to the current group
                    currentGroup.push_back(hash_b[i].second);
                } else {
                    // If the current element is different, start a new group
                    if(currentGroup.size() > 1){
                        m_twins.push_back(currentGroup);
                    }
                    currentGroup.clear();
                    currentGroup.push_back(hash_b[i].second);
                }
            }
            // Add the last group
            if(currentGroup.size() > 1){
                m_twins.push_back(currentGroup);
            }

            // sort and make unique
            for (auto &twins: m_twins) {
                std::sort(twins.begin(), twins.end());
                twins.erase(std::unique(twins.begin(), twins.end()), twins.end());
            }

            for (auto &twins: m_twins) {
                ASSERT(std::is_sorted(twins.begin(), twins.end()));
                ASSERT(no_duplicates(twins));
                ASSERT(twins.size() >= 2);
            }
            ASSERT(no_duplicates(m_twins));

            for (u32 i = 0; i < (u32) m_twins.size(); ++i) {
                for (u32 j = i + 1; j < (u32) m_twins.size(); ++j) {
                    for (u32 k = 0; k < (u32) m_twins[i].size(); ++k) {
                        ASSERT(!exists(m_twins[j], m_twins[i][k]));
                    }
                }
            }

        }

        /**
         * Removes all twins from the graph and creates the resulting graph.
         *
         * @param g The graph to be freed.
         * @return The graph free of twins.
         */
        inline Graph reduce_twins(const Graph &g) {
            u32 n_A = g.n_A;
            u32 n_B = g.n_B;

            // adjust the new number of vertices in B
            for (const auto &twins: m_twins) {
                n_B = (n_B + 1) - (u32) twins.size();
            }

            Graph new_g(n_A, n_B);

            u32 new_vertex_b = 0;
            for (u32 vertex_b = 0; vertex_b < g.n_B; ++vertex_b) {
                // process each vertex of B

                // check if vertex belongs to a twin, and if it is the smallest vertex
                bool vertex_found = false;
                bool vertex_smallest = false;
                u32 twins_size = 1;
                for (auto &twins: m_twins) {

                    ASSERT(std::is_sorted(twins.begin(), twins.end()));
                    ASSERT(no_duplicates(twins));
                    ASSERT(twins.size() >= 2);

                    if (exists(twins, vertex_b)) {
                        vertex_found = true;
                        vertex_smallest = (vertex_b == min(twins));
                        twins_size = (u32) twins.size();
                    }
                }

                if ((vertex_found && vertex_smallest) || !vertex_found) {
                    // process this vertex
                    m_tt_twins.add_B(vertex_b, new_vertex_b);

                    for (Edge vertex_a: g.adj_list[vertex_b]) {
                        new_g.add_edge(vertex_a.vertex, new_vertex_b, vertex_a.weight * twins_size);
                    }
                    new_vertex_b += 1;
                }
            }

            new_g.finalize();
            return new_g;
        }

        /**
         * Back propagates the solution for the twin reduction.
         *
         * @param sol The solution of the twin reduced graph.
         * @return The solution before the twin reduced graph.
         */
        inline AlignedVector<u32> back_propagate_twins(const AlignedVector<u32> &sol) {
            AlignedVector<u32> new_sol;

            for (u32 vertex: sol) {
                u32 old_vertex = m_tt_twins.get_B_old(vertex);

                // check if vertex belongs to a twin, and if it is the smallest vertex
                bool is_twin = false;
                for (auto &twins: m_twins) {

                    ASSERT(std::is_sorted(twins.begin(), twins.end()));
                    ASSERT(no_duplicates(twins));
                    ASSERT(twins.size() >= 2);

                    if (exists(twins, old_vertex)) {
                        // insert the whole twin
                        is_twin = true;
                        for (u32 twin: twins) { new_sol.push_back(twin); }
                    }
                }

                if (!is_twin) {
                    new_sol.push_back(old_vertex);
                }
            }

            return new_sol;
        }
    };

}

#endif //PACE2024EXACT_TWIN_REDUCER_H
