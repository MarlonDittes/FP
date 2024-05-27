#ifndef PACE2024EXACT_ONE_VERTEX_PARTITIONER_H
#define PACE2024EXACT_ONE_VERTEX_PARTITIONER_H

#include "definitions.h"
#include "macros.h"
#include "misc.h"
#include "graph_hen.h"
#include "translation_table.h"

namespace CrossGuard {

    /**
     * Assume the scenario:
     * AAAAAAAAAAAAAAAA C   G DDDDDDDDDDDDDDDDDD
     * xxxxxxxxxxxxxxxx  \ /  xxxxxxxxxxxxxxxxxx
     * BBBBBBBBBBBBBBBB   E   FFFFFFFFFFFFFFFFFF
     *
     * AC and B form one component, and GD and F form one component. Both
     * components are only connected via the vertex E. Solving (AC,B) and (GD,F)
     * is simpler than solving the whole graph.
     */
    class OneVertexPartitioner {
    public:
        const Graph &m_graph;

        AlignedVector<int> m_component_id;

        AlignedVector<u32> m_left_A;
        AlignedVector<u32> m_left_B;
        TranslationTable m_left_tt;

        AlignedVector<u32> m_right_A;
        AlignedVector<u32> m_right_B;
        TranslationTable m_right_tt;

        u32 connecting_vertex = std::numeric_limits<u32>::max();
        u32 n_comp = 1;

        explicit OneVertexPartitioner(const Graph &g) : m_graph(g) {
            m_component_id.resize(m_graph.n_B);
        }

        /**
         * Finds all components and creates sub-graphs with translation tables.
         */
        inline void find_components() {
            identify_components();
        }

        /**
         * Returns all components of the graph.
         *
         * @return Vector with sub-graphs.
         */
        inline AlignedVector<Graph> get_components() {
            AlignedVector<Graph> new_graphs;

            // if no components, simply return the original graph
            bool one_graph_empty = m_left_B.empty() || m_right_B.empty();
            if (connecting_vertex == std::numeric_limits<u32>::max() || one_graph_empty) {
                new_graphs.push_back(m_graph);
                return new_graphs;
            }

            {
                // left graph
                u32 n_A = m_left_A.size();
                u32 n_B = m_left_B.size();
                Graph g_left(n_A, n_B);

                // translate
                for (u32 i = 0; i < m_left_A.size(); ++i) {
                    m_left_tt.add_A(m_left_A[i], i);
                }
                for (u32 i = 0; i < m_left_B.size(); ++i) {
                    m_left_tt.add_B(m_left_B[i], i);
                }

                // add edges
                for (u32 old_b: m_left_B) {
                    for (Edge a: m_graph.adj_list[old_b]) {
                        u32 new_b = m_left_tt.get_B_new(old_b);
                        u32 new_a = m_left_tt.get_A_new(a.vertex);
                        u32 weight = a.weight;
                        g_left.add_edge(new_a, new_b, weight);
                    }
                }

                g_left.finalize();
                new_graphs.push_back(g_left);
            }
            {
                // right graph
                u32 n_A = m_right_A.size();
                u32 n_B = m_right_B.size();
                Graph g_right(n_A, n_B);

                // translate
                for (u32 i = 0; i < m_right_A.size(); ++i) {
                    m_right_tt.add_A(m_right_A[i], i);
                }
                for (u32 i = 0; i < m_right_B.size(); ++i) {
                    m_right_tt.add_B(m_right_B[i], i);
                }

                // add edges
                for (u32 old_b: m_right_B) {
                    for (Edge a: m_graph.adj_list[old_b]) {
                        u32 new_b = m_right_tt.get_B_new(old_b);
                        u32 new_a = m_right_tt.get_A_new(a.vertex);
                        u32 weight = a.weight;
                        g_right.add_edge(new_a, new_b, weight);
                    }
                }

                g_right.finalize();
                new_graphs.push_back(g_right);
            }

            return new_graphs;
        }

        /**
         * Back propagates all solutions of the sub graphs.
         *
         * @param sols The solutions of the sub graphs.
         * @param order The order of the solutions.
         * @return The solution of the original graph.
         */
        inline AlignedVector<u32> back_propagate(AlignedVector<AlignedVector<u32>> &sols) const {
            AlignedVector<u32> new_sol;

            if (sols.size() == 1) {
                new_sol = sols[0];
                return new_sol;
            }

            ASSERT(sols.size() == 2);
            for (u32 new_b: sols[0]) {
                u32 old_b = m_left_tt.get_B_old(new_b);
                new_sol.push_back(old_b);
            }
            new_sol.push_back(connecting_vertex);
            for (u32 new_b: sols[1]) {
                u32 old_b = m_right_tt.get_B_old(new_b);
                new_sol.push_back(old_b);
            }

            return new_sol;
        }

    private:
        void identify_components() {
            AlignedVector<u32> queue;
            for (u32 i = 0; i < m_graph.n_B; ++i) {
                // std::cout << "i = " << i << std::endl;
                // print(m_graph.adj_list[i]);
                if (can_be_connecting_vertex(i)){
                    // determine if this is a special vertex
                    ASSERT(m_graph.adj_list[i][0].vertex < m_graph.adj_list[i][1].vertex);
                    queue.clear();
                    std::fill(m_component_id.begin(), m_component_id.end(), -1);

                    // process the "left" part
                    for (u32 a = 0; a <= m_graph.adj_list[i][0].vertex; ++a) {
                        for (Edge b: m_graph.adj_list_from_A[a]) {
                            m_component_id[b.vertex] = 0;
                        }
                    }

                    // std::cout << "i = " << i << std::endl;
                    // print(m_component_id);

                    // check the right part
                    bool right_part_touched = false;
                    for (u32 a = m_graph.adj_list[i][1].vertex; a < m_graph.n_A; ++a) {
                        for (Edge b: m_graph.adj_list_from_A[a]) {
                            if (b.vertex != i && m_component_id[b.vertex] == 0) {
                                right_part_touched = true;
                                a = std::numeric_limits<u32>::max() - 1;
                                break;
                            }
                        }
                    }

                    // check if the right part has been touched
                    if (right_part_touched) {
                        // right part has been touched, this is not a special vertex
                        continue;
                    } else {
                        // the right vertex in A has not been touched, this is a special vertex.
                        connecting_vertex = i;
                        n_comp = 2;
                        // std::cout << "yes" << std::endl;

                        for (u32 j = 0; j < m_graph.n_B; ++j) {
                            if (j != connecting_vertex) {
                                // do not add the connecting vertex
                                if (m_component_id[j] == 0) {
                                    m_left_B.push_back(j);
                                } else {
                                    m_right_B.push_back(j);
                                }
                            }
                        }

                        for (u32 b: m_left_B) {
                            for (Edge a: m_graph.adj_list[b]) {
                                m_left_A.push_back(a.vertex);
                            }
                        }
                        std::sort(m_left_A.begin(), m_left_A.end());
                        m_left_A.erase(std::unique(m_left_A.begin(), m_left_A.end()), m_left_A.end());

                        for (u32 b: m_right_B) {
                            for (Edge a: m_graph.adj_list[b]) {
                                m_right_A.push_back(a.vertex);
                            }
                        }
                        std::sort(m_right_A.begin(), m_right_A.end());
                        m_right_A.erase(std::unique(m_right_A.begin(), m_right_A.end()), m_right_A.end());

                        break;
                    }
                }
            }
        }

        bool can_be_connecting_vertex(u32 b){
            if(m_graph.adj_list[b].size() < 2){
                return false;
            }

            // determine if all vertices of A are connected only to b
            for(u32 i = 2; i < m_graph.adj_list[b].size(); ++i){
                Edge a = m_graph.adj_list[b][i - 1];
                if(!(m_graph.adj_list_from_A[a.vertex].size() == 1 && m_graph.adj_list_from_A[a.vertex][0].vertex == b)){
                    return false;
                }
            }

            // determine if all vertices of A are in one line
            for(u32 i = 0; i < m_graph.adj_list[b].size() - 1; ++i){
                if(m_graph.adj_list[b][i].vertex + 1 != m_graph.adj_list[b][i + 1].vertex){
                    return false;
                }
            }

            return true;
        }
    };

}

#endif //PACE2024EXACT_ONE_VERTEX_PARTITIONER_H
