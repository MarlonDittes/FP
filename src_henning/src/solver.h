#ifndef PACE2024EXACT_SOLVER_H
#define PACE2024EXACT_SOLVER_H

#include "graph_hen.h"
#include "partitioner.h"
#include "exhaustive_solver.h"
#include "twin_reducer.h"
#include "useless_reducer.h"
#include "front_back_reducer.h"
#include "one_vertex_partitioner.h"
#include "one_vertex_partition_manager.h"

namespace CrossGuard {

    // TODO: How to include Domination Reduction
    // TODO: Is there a way to handle vertices b with only one neighbor?


/**
 * Class to solve the problem.
 */
    class Solver {
    private:
        const Graph &m_graph;
        AlignedVector<u32> m_solution;

        // timings
        std::chrono::steady_clock::time_point sp;
        std::chrono::steady_clock::time_point ep;

    public:
        /**
         * Default constructor.
         *
         * @param g The graph.
         */
        explicit Solver(const Graph &g) : m_graph(g) {}

        /**
         * Solves the problem.
         */
        inline void solve(bool verbose=false) {
            sp = std::chrono::steady_clock::now();
            ASSERT(m_graph.is_finalized);

            // find components of the graph
            Partitioner partitioner(m_graph);
            partitioner.find_components();

            // determine the component order
            Graph partition_g = partitioner.get_component_graph();
            ASSERT(partition_g.is_finalized);
            ExhaustiveSolver component_solver(partition_g);
            component_solver.solve();
            AlignedVector<u32> component_order = component_solver.get_solution();

            AlignedVector<Graph> partitions = partitioner.get_components();

            if(verbose){
                AlignedVector<u32> sizes;
                for(auto &g_part : partitions){sizes.push_back(g_part.n_B);}
                std::cout << "\tNumber of components: " << component_order.size() << " ";
                print(sizes);
            }

            // solve each component
            AlignedVector<AlignedVector<u32>> solutions;
            for (auto &g: partitions) {
                // solve sub-graph
                ASSERT(g.is_finalized);
                // g.print();

                // further partition the graph
                OneVertexPartitionManager one_vertex_partition_manager(g);
                AlignedVector<Graph> new_partitions = one_vertex_partition_manager.get_components();

                if(verbose){
                    AlignedVector<u32> sizes;
                    for(auto &g_part : new_partitions){sizes.push_back(g_part.n_B);}
                    std::cout << "\t\tFurther divided: " << new_partitions.size() << " ";
                    print(sizes);
                }

                AlignedVector<AlignedVector<u32>> one_vertex_solutions;
                for(auto &g_part : new_partitions) {
                    ASSERT(g_part.is_finalized);
                    // g_part.print();

                    UselessReducer useless_reducer(g_part);
                    Graph useless_g = useless_reducer.reduce();
                    ASSERT(useless_g.is_finalized);

                    FrontBackReducer fb_reducer(useless_g);
                    Graph fb_g = fb_reducer.reduce();
                    ASSERT(fb_g.is_finalized);

                    TwinReducer reducer(fb_g);
                    Graph reduced_g = reducer.reduce();
                    ASSERT(reduced_g.is_finalized);

                    if (verbose) {
                        std::cout << "\t\t\tComponent reduced from " << g_part.n_B << " to " << useless_g.n_B << " to " << fb_g.n_B << " to " << reduced_g.n_B << std::flush;
                    }

                    ExhaustiveSolver s(reduced_g);
                    AlignedVector<u32> median_vector = reduced_g.get_median_solution();
                    s.set_initial_solution(median_vector);
                    s.solve();
                    AlignedVector<u32> temp = s.get_solution();

                    if (verbose) { std::cout << " -- Component solved" << std::endl; }

                    temp = reducer.back_propagate(temp);
                    temp = fb_reducer.back_propagate(temp);
                    temp = useless_reducer.back_propagate(temp);
                    one_vertex_solutions.push_back(temp);
                }
                AlignedVector<u32> temp2 = one_vertex_partition_manager.back_propagate(one_vertex_solutions);
                solutions.push_back(temp2);
            }
            m_solution = partitioner.back_propagate(solutions, component_order);

            ep = std::chrono::steady_clock::now();
        }

        /**
         * Returns the solution.
         *
         * @return The solution.
         */
        inline AlignedVector<unsigned int> get_solution() {
            AlignedVector<u32> sol(m_solution);
            return sol;
        }

        /**
         * Returns the permutation vector. All entries are in the range
         * [n_A + 1, ..., n_A + n_B].
         *
         * @return Permutation of B.
         */
        inline AlignedVector<unsigned int> get_shifted_solution() const {
            AlignedVector<u32> v(m_solution.size());
            std::copy(m_solution.begin(), m_solution.end(), v.begin());
            for (auto &x: v) {
                x += m_graph.n_A + 1;
            }
            return v;
        }

        inline double get_time() {
            return get_seconds(sp, ep);
        }
    };

}

#endif //PACE2024EXACT_SOLVER_H
