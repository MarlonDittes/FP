#ifndef PACE2024EXACT_FRONT_BACK_REDUCER_H
#define PACE2024EXACT_FRONT_BACK_REDUCER_H

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
     * Uses Front Back Reduction to fix vertices either at the front or the back.
     */
    class FrontBackReducer {
    private:
        const Graph &m_graph;

        AlignedVector<u32> m_front;
        AlignedVector<u32> m_back;

        TranslationTable m_tt;

    public:
        /**
         * Default constructor.
         *
         * @param g The graph.
         */
        explicit FrontBackReducer(const Graph &g) : m_graph(g) {}

        /**
         * Reduces the graph and returns a new graph.
         *
         * @return The newly reduced graph.
         */
        inline Graph reduce() {
            ASSERT(m_graph.is_finalized);
            Graph g = m_graph;

            find_front_back_reductions();
            g = reduce_front_back_reductions();

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
            AlignedVector<u32> new_sol;
            new_sol.reserve(m_graph.n_B);

            for(u32 i : m_front){
                new_sol.push_back(i);
            }
            for(u32 i : sol){
                u32 old_i = m_tt.get_B_old(i);
                new_sol.push_back(old_i);
            }
            for(u32 i : m_back){
                new_sol.push_back(i);
            }

            return new_sol;
        }

    private:

        void find_front_back_reductions(){
            // determine for each vertex, whether all its connections are at the front of A
            for(u32 i = 0; i < m_graph.n_B; ++i){

                if(m_graph.adj_list[i].size() > 2){
                    continue;
                }

                bool all_front = true;
                ASSERT(std::is_sorted(m_graph.adj_list[i].begin(), m_graph.adj_list[i].end()));
                for(u32 j = 0; j < m_graph.adj_list[i].size(); ++j){
                    all_front &= (j == m_graph.adj_list[i][j].vertex);
                }

                if(all_front){
                    m_front.push_back(i);
                }
            }

            // determine for each vertex, whether all its connections are at the back of A
            for(u32 i = 0; i < m_graph.n_B; ++i){

                if(m_graph.adj_list[i].size() > 2){
                    continue;
                }

                bool all_back = true;
                ASSERT(std::is_sorted(m_graph.adj_list[i].begin(), m_graph.adj_list[i].end()));
                for(u32 j = 0; j < m_graph.adj_list[i].size(); ++j){
                    all_back &= (m_graph.n_A - 1 - j == m_graph.adj_list[i][m_graph.adj_list[i].size() - 1 - j].vertex);
                }

                if(all_back && !exists(m_front, i)){
                    m_back.push_back(i);
                }
            }

            // sorts both ends by neighborhood size
            std::sort(m_front.begin(), m_front.end(), [&](u32 i, u32 j){return m_graph.adj_list[i].size() < m_graph.adj_list[j].size();});
            std::sort(m_back.begin(), m_back.end(), [&](u32 i, u32 j){return m_graph.adj_list[i].size() > m_graph.adj_list[j].size();});
            // std::cout << "front and back" << std::endl;
            // print(m_front);
            // print(m_back);
        }

        Graph reduce_front_back_reductions(){
            // in case all vertices can be reduced, leave one trivial vertex
            if(m_front.size() + m_back.size() == m_graph.n_B){
                if(!m_front.empty()){
                    m_front.pop_back();
                } else {
                    // shift one to the front
                    for(u32 i = 0; i < m_back.size() - 1; ++i){
                        m_back[i] = m_back[i + 1];
                    }
                    m_back.pop_back();
                }
            }

            // std::cout << "front and back" << std::endl;
            // print(m_front);
            // print(m_back);

            u32 n_A = m_graph.n_A;
            u32 n_B = m_graph.n_B - m_front.size() - m_back.size();
            Graph new_g(n_A, n_B);

            u32 new_i = 0;
            for(u32 i = 0; i < m_graph.n_B; ++i){
                if (!exists(m_front, i) && !exists(m_back, i)){
                    // process i
                    for(Edge vertex_a : m_graph.adj_list[i]){
                        new_g.add_edge(vertex_a.vertex, new_i, vertex_a.weight);
                    }
                    m_tt.add_B(i, new_i);
                    new_i += 1;
                }
            }

            new_g.finalize();
            return new_g;
        }
    };

}

#endif //PACE2024EXACT_FRONT_BACK_REDUCER_H
