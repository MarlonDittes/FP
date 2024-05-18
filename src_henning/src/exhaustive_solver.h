#ifndef PACE2024EXACT_EXHAUSTIVESOLVER_H
#define PACE2024EXACT_EXHAUSTIVESOLVER_H

#include <vector>
#include <cstdint>
#include <numeric>
#include <limits>

#include "definitions.h"
#include "macros.h"
#include "misc.h"
#include "graph_hen.h"
#include "cross_matrix.h"
#include "candidate_manager.h"
#include "perfect_pattern.h"

namespace CrossGuard {

    struct DepthInfo {
        s32 counter;
        u32 curr_cuts; // number of cuts between placed vertices
        u32 min_intra_cuts; // min number of cuts between non-placed vertices
        u32 min_inter_cuts; // min number of cuts between placed and non-placed vertices
    };

    /**
     * Solver for the one-sided crossing minimization problem (OCM).
     */
    class ExhaustiveSolver {
    private:
        const Graph &m_graph;
        const u32 n;

        AlignedVector<CandidateManager> m_candidate_order;
        AlignedVector<DepthInfo> m_depth_info;

        // Minimum Cross Matrix
        CrossMatrix m_cross_matrix;
        PerfectPattern m_perfect_pattern;

        // holds the best found permutation
        AlignedVector<u32> m_solution;
        u32 m_solution_n_cuts;

        // min number of cuts possible
        u32 m_min_n_cuts;

    public:
        /**
         * Default constructor.
         *
         * @param g The graph to optimize.
         */
        explicit ExhaustiveSolver(Graph &g) : m_graph(g),
                                              n(g.n_B),
                                              m_cross_matrix(g.n_B),
                                              m_perfect_pattern(g.n_B) {
            m_candidate_order.reserve(n);
            for (u32 i = 0; i < n; ++i) {
                m_candidate_order.emplace_back(1 + n - i);
            }

            m_depth_info.resize(n);

            m_solution.resize(n); // reserve space
            m_solution_n_cuts = std::numeric_limits<u32>::max();

            // initialize the Crossing Matrix
            m_min_n_cuts = std::numeric_limits<u32>::max();
            initialize_CrossMatrix();
        }

        /**
         * Sets an initial solution.
         *
         * @param init_sol The initial solution.
         */
        inline void set_initial_solution(AlignedVector<unsigned int> &init_sol) {
            check_new_solution(init_sol);
        }

        /**
         * Determines the m_permutation, with the least number of cuts.
         */
        inline void solve() {
            // most simple solution if required
            if (m_solution_n_cuts == std::numeric_limits<u32>::max()) {
                std::iota(m_solution.begin(), m_solution.end(), 0);
                check_new_solution(m_solution);
            }

            // exhaustive search
            iterative_search();
        }

        /**
         * Returns the permutation vector. All entries are in the range
         * [0, ..., n - 1].
         *
         * @return Permutation of B.
         */
        inline AlignedVector<unsigned int> get_solution() const {
            AlignedVector<unsigned int> v(m_solution.size());
            std::copy(m_solution.begin(), m_solution.end(), v.begin());
            return v;
        }

    private:
        /**
         * Determines is a new solution has been found.
         *
         * @param sol The possible solution.
         */
        inline bool check_new_solution(AlignedVector<u32> &sol) {
            u32 n_cuts = count_cuts(sol, n);
            if (n_cuts < m_solution_n_cuts) {
                std::copy(sol.begin(), sol.end(), m_solution.begin());
                m_solution_n_cuts = n_cuts;
                return true;
            }
            return false;
        }

        /**
         * Counts the number of cuts, based on the current permutation.
         *
         * @return Number of cuts.
         */
        inline u32 count_cuts(AlignedVector<u32> &perm, u32 size) {
            u32 n_cuts = 0;
            for (u32 i = 0; i < size; ++i) {
                for (u32 j = i + 1; j < size; ++j) {
                    n_cuts += m_cross_matrix.matrix[perm[i] * n + perm[j]];
                }
            }
            return n_cuts;
        }

        /**
         * Initialize the Cross-Matrix.
         */
        inline void initialize_CrossMatrix() {
            ASSERT(m_graph.is_finalized);
            m_min_n_cuts = 0;

            for (u32 i = 0; i < n; ++i) {
                for (u32 j = i + 1; j < n; ++j) {
                    u32 c1 = 0;
                    u32 c2 = 0;
                    // loop through the edges
                    for (size_t k = 0; k < m_graph.adj_list[i].size(); ++k) {
                        for (size_t l = 0; l < m_graph.adj_list[j].size(); ++l) {
                            Edge a1 = m_graph.adj_list[i][k];
                            Edge a2 = m_graph.adj_list[j][l];

                            u32 cuts = a1.weight * a2.weight;

                            c1 += (a2.vertex < a1.vertex) * cuts;
                            c2 += (a1.vertex < a2.vertex) * cuts;
                        }
                    }
                    m_cross_matrix.matrix[i * n + j] = c1;
                    m_cross_matrix.matrix[j * n + i] = c2;

                    m_min_n_cuts += std::min(c1, c2);

                    if (c1 == 0 && c2 > 0) { m_perfect_pattern.add(i, j); }
                    if (c2 == 0 && c1 > 0) { m_perfect_pattern.add(j, i); }
                }
            }
        }

        /**
         * Iteratively searches for a solution.
         */
        inline void iterative_search() {
            execute_DCO_depth_0();
            m_depth_info[0] = {-1, 0, calculate_intraMinCuts_depth_0(), 0};

            s32 depth = 0;

            while (depth >= 0) {

                // move further
                m_depth_info[depth].counter += 1;

                // check if at the end
                if (m_depth_info[depth].counter == ((s32) n - depth)) {
                    depth -= 1;
                    continue;
                }

                // check if we arrived at the bottom
                if (depth == (s32) n - 1) {
                    u32 last_vertex_gain = m_candidate_order[depth].candidates[m_depth_info[depth].counter].gain;

                    if (m_depth_info[depth].curr_cuts + last_vertex_gain < m_solution_n_cuts) {
                        // transfer solution
                        for (u32 i = 0; i < n; ++i) {
                            m_solution[i] = m_candidate_order[i].candidates[m_depth_info[i].counter].c;
                        }
                        m_solution_n_cuts = m_depth_info[depth].curr_cuts + last_vertex_gain;

                        if (m_solution_n_cuts == m_min_n_cuts) {
                            // there is no better solution
                            return;
                        }

                        // we updated the best solution, check which depth to jump to
                        while (depth >= 0 && m_depth_info[depth].curr_cuts + m_depth_info[depth].min_intra_cuts + m_depth_info[depth].min_inter_cuts >= m_solution_n_cuts) {
                            depth -= 1;
                        }
                        depth += 1;
                    }

                    // go up
                    depth -= 1;
                    continue;
                }

                // move down and reset counter
                depth += 1;
                // reset the counter
                m_depth_info[depth].counter = -1;

                bool abort = execute_DCO(depth);
                if (abort) {
                    // go up
                    depth -= 1;
                }
            }
        }

        /**
         * Executes Dynamic-Candidate-Ordering at depth 0.
         */
        inline void execute_DCO_depth_0() {
            m_candidate_order[0].initialize_median_vec(n);
            for (u32 i = 0; i < n; ++i) {
                m_candidate_order[0].candidates[i] = {i, 0};
                m_candidate_order[0].median_candidate[i] = {i, m_graph.medians[i]};
                m_candidate_order[0].max_median = std::max(m_candidate_order[0].max_median, m_graph.medians[i]);
            }
            m_candidate_order[0].sort(n);
        };

        /**
         * Executes Dynamic-Candidate-Ordering.
         *
         * @param depth Current depth of the algorithm.
         */
        inline bool execute_DCO(int depth) {
            // get last added vertex
            Candidate last_added_vertex = m_candidate_order[depth - 1].candidates[m_depth_info[depth - 1].counter];

            // determine the number of cuts between placed and still to be placed vertices
            u32 inter_cuts = 0;

            // count the number of cuts that get removed by adding the last vertex
            u32 n_intra_cuts_removed = 0;

            u32 idx = 0;
            for (u32 i = 0; i < n - (depth - 1); ++i) {
                // curr vertex
                Candidate vertex = m_candidate_order[depth - 1].candidates[i];

                u32 vertex_not_used = (1 - (vertex.c == last_added_vertex.c)); // 0 if vertex is used, 1 else

                // check if a pattern is violated
                if (m_perfect_pattern.matches(vertex.c, last_added_vertex.c)) {
                    return true;
                }

                // update inter cuts
                u32 extra_gain = m_cross_matrix.matrix[last_added_vertex.c * n + vertex.c];
                inter_cuts += vertex_not_used * (vertex.gain + extra_gain);

                // update intra cuts
                n_intra_cuts_removed += vertex_not_used * std::min(extra_gain, m_cross_matrix.matrix[vertex.c * n + last_added_vertex.c]);

                // place vertex
                m_candidate_order[depth].candidates[idx] = vertex;
                m_candidate_order[depth].candidates[idx].gain += extra_gain;

                idx += vertex_not_used;
            }

            // update cuts
            m_depth_info[depth].curr_cuts = m_depth_info[depth - 1].curr_cuts + last_added_vertex.gain;
            m_depth_info[depth].min_intra_cuts = m_depth_info[depth - 1].min_intra_cuts - n_intra_cuts_removed;
            m_depth_info[depth].min_inter_cuts = inter_cuts;

            // check if we can abort based on the number of cuts remaining
            if (m_depth_info[depth].curr_cuts + m_depth_info[depth].min_intra_cuts + m_depth_info[depth].min_inter_cuts >= m_solution_n_cuts) {
                return true;
            }

            return false;
        }

        /**
         * Gets the minimum amount of cuts between the candidates still to be
         * placed.
         *
         * @return Minimum number of cuts.
         */
        inline u32 calculate_intraMinCuts_depth_0() {
            u32 n_cuts = 0;
            for (u32 i = 0; i < n; ++i) {
                for (u32 j = i + 1; j < n; ++j) {
                    u32 c1 = m_candidate_order[0].candidates[i].c;
                    u32 c2 = m_candidate_order[0].candidates[j].c;

                    u32 cut1 = m_cross_matrix.matrix[c1 * n + c2];
                    u32 cut2 = m_cross_matrix.matrix[c2 * n + c1];

                    n_cuts += std::min(cut1, cut2);
                }
            }
            return n_cuts;
        }
    };

}


#endif //PACE2024EXACT_EXHAUSTIVESOLVER_H
