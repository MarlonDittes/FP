#ifndef PACE2024EXACT_SOLVER_BF_H
#define PACE2024EXACT_SOLVER_BF_H

#define BF_DEBUG 0

#include <vector>
#include <cstdint>
#include <numeric>
#include <limits>

#include "definitions.h"
#include "macros.h"
#include "misc.h"
#include "graph_hen.h"
#include "cross_matrix.h"

namespace CrossGuard {

/**
 * Brute-Force Solver.
 */
    class Solver_BF {
    private:
        const Graph &m_graph;

        AlignedVector<u32> permutation; // current linear order
        AlignedVector<bool> is_used; // O(1) access if a vertex is already used
        u32 curr_size;

        AlignedVector<u32> solution;
        u32 solution_n_cuts;

        // Minimum Cross Matrix
        CrossMatrix m_cross_matrix;

        // debug vars
#if BF_DEBUG
        AlignedVector<AlignedVector<u32>> all_permutations;
#endif

    public:

        /**
         * Default constructor.
         *
         * @param g The m_graph to optimize.
         */
        explicit Solver_BF(Graph &g) : m_graph(g), m_cross_matrix(m_graph.n_B) {
            permutation.resize(m_graph.n_B);
            is_used.resize(m_graph.n_B, false);
            curr_size = 0;

            solution.resize(m_graph.n_B);
            solution_n_cuts = std::numeric_limits<u32>::max();
        }

        /**
         * Determines the m_permutation, with the least number of cuts.
         */
        inline void solve() {
            initialize_CrossMatrix();

            recursive_solve();

#if BF_DEBUG
            // check that n! permutations have been found
            if ((u32) all_permutations.size() != fac(m_graph.n_B)) {
                std::cout << "WARNING: Not all permutations have been found! Only " << all_permutations.size() << " of " << fac(m_graph.n_B) << " have been found!" << std::endl;
            }

            // check that each entry has exactly m_n_B entries
            for (auto &vec: all_permutations) {
                if ((u32) vec.size() != m_graph.n_B) {
                    std::cout << "WARNING: Permutation does not contain " << m_graph.n_B << " elements, but " << vec.size() << "!" << std::endl;
                    std::cout << "Permutation: ";
                    print(vec);
                }
            }

            // check that no entry has duplicates
            for (auto &vec: all_permutations) {
                if (!no_duplicates(vec)) {
                    std::cout << "WARNING: Permutation contains duplicates!" << std::endl;
                    std::cout << "Permutation: ";
                    print(vec);
                }
            }

            // check that same permutation has not been found multiple times
            if(!no_duplicates(all_permutations)){
                std::cout << "WARNING: Permutation has been considered multiple times!" << std::endl;
            }
#endif
        }

        /**
         * Returns the m_permutation vector. All entries are in the range
         * [0, ..., n_B - 1].
         *
         * @return Permutation of B.
         */
        inline AlignedVector<unsigned int> get_solution() const {
            AlignedVector<unsigned int> v(solution.size());
            for(size_t i = 0; i < solution.size(); ++i){
                v[i] = solution[i];
            }
            return v;
        }

        /**
         * Returns the m_permutation vector. All entries are in the range
         * [n_A + 1, ..., n_A + n_B].
         *
         * @return Permutation of B.
         */
        inline AlignedVector<unsigned int> get_shifted_solution() const {
            AlignedVector<unsigned int> v(solution.size());
            for(size_t i = 0; i < solution.size(); ++i){
                v[i] = solution[i] + m_graph.n_A + 1;
            }
            return v;
        }


    private:
        /**
         * Counts the number of cuts, based on the current m_permutation.
         *
         * @return Number of cuts.
         */
        inline u32 count_cuts() {
            u32 n_cuts = 0;
            for (u32 i = 0; i < curr_size; ++i) {
                for (u32 j = i + 1; j < curr_size; ++j) {
                    n_cuts += m_cross_matrix.matrix[permutation[i] * m_graph.n_B + permutation[j]];
                }
            }
            return n_cuts;
        }

        /**
         * Initialize the Cross-Matrix.
         */
        inline void initialize_CrossMatrix() {
            ASSERT(m_graph.is_finalized);

            for (u32 i = 0; i < m_graph.n_B; ++i) {
                for (u32 j = i+1; j < m_graph.n_B; ++j) {
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
                    m_cross_matrix.matrix[i * m_graph.n_B + j] = c1;
                    m_cross_matrix.matrix[j * m_graph.n_B + i] = c2;
                }
            }
        }


        /**
         * Recursively searches the permutation tree.
         */
        inline void recursive_solve() {
            if (curr_size == m_graph.n_B) {
#if BF_DEBUG
                // add found m_permutation
                all_permutations.push_back(permutation);
#endif
                // we have a permutation, check the number of cuts
                u32 n_cuts = count_cuts();
                if (n_cuts < solution_n_cuts) {
                    std::copy(permutation.begin(), permutation.end(), solution.begin());
                    solution_n_cuts = n_cuts;
                }
                return;
            }

            for (u32 i = 0; i < m_graph.n_B; ++i) {
                if (!is_used[i]) {
                    // vertex i is not used, so append it to the m_permutation
                    permutation[curr_size] = i;
                    is_used[i] = true;
                    curr_size += 1;

                    recursive_solve();

                    is_used[i] = false;
                    curr_size -= 1;
                }
            }
        }
    };

}

#endif //PACE2024EXACT_SOLVER_BF_H
