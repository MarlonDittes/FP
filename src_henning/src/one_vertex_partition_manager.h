#ifndef PACE2024EXACT_ONE_VERTEX_PARTITION_MANAGER_H
#define PACE2024EXACT_ONE_VERTEX_PARTITION_MANAGER_H

#include "definitions.h"
#include "macros.h"
#include "misc.h"
#include "graph_hen.h"
#include "translation_table.h"
#include "one_vertex_partitioner.h"

namespace CrossGuard {

    /**
     * Recursively applies one vertex partitioning.
     */
    class OneVertexPartitionManager{
    public:
        const Graph &m_graph;

        AlignedVector<Graph> m_graphs;

        AlignedVector<u32> m_connecting_vertices;
        AlignedVector<TranslationTable> m_tts;

        explicit OneVertexPartitionManager(const Graph &g) : m_graph(g) {}

        AlignedVector<Graph> get_components(){
            TranslationTable init_tt;
            for(u32 i = 0; i < m_graph.n_A; ++i){
                init_tt.add_A(i, i);
            }
            for(u32 i = 0; i < m_graph.n_B; ++i){
                init_tt.add_B(i, i);
            }

            process(m_graph, init_tt);

            return m_graphs;
        }

        AlignedVector<u32> back_propagate(AlignedVector<AlignedVector<u32>> &solutions){
            AlignedVector<u32> new_sol;

            ASSERT(solutions.size() == m_graphs.size());
            ASSERT(solutions.size() == m_tts.size());
            ASSERT(m_connecting_vertices.size() == solutions.size() - 1);

            for(u32 i = 0; i < m_connecting_vertices.size(); ++i){
                for(u32 x : solutions[i]){
                    new_sol.push_back(m_tts[i].get_B_old(x));
                }
                new_sol.push_back(m_connecting_vertices[i]);
            }
            for(u32 x : solutions.back()){
                new_sol.push_back(m_tts.back().get_B_old(x));
            }

            return new_sol;
        }

    private:
        void process(const Graph &g, TranslationTable &tt){
            OneVertexPartitioner partitioner(g);
            partitioner.find_components();
            AlignedVector<Graph> components = partitioner.get_components();

            if(components.size() == 1){
                // end of recursion
                m_graphs.push_back(g);
                m_tts.push_back(tt);

                return;
            }
            ASSERT(components.size() == 2);

            // create left translation table
            TranslationTable left_tt;
            for(u32 a = 0; a < components[0].n_A; ++a){
                left_tt.add_A(tt.get_A_old(partitioner.m_left_tt.get_A_old(a)), a);
            }
            for(u32 b = 0; b < components[0].n_B; ++b){
                left_tt.add_B(tt.get_B_old(partitioner.m_left_tt.get_B_old(b)), b);
            }
            process(components[0], left_tt);

            // create right translation table
            TranslationTable right_tt;
            for(u32 a = 0; a < components[1].n_A; ++a){
                right_tt.add_A(tt.get_A_old(partitioner.m_right_tt.get_A_old(a)), a);
            }
            // add the connecting vertex
            m_connecting_vertices.push_back(tt.get_B_old(partitioner.connecting_vertex));
            for(u32 b = 0; b < components[1].n_B; ++b){
                right_tt.add_B(tt.get_B_old(partitioner.m_right_tt.get_B_old(b)), b);
            }
            process(components[1], right_tt);
        }
    };

}

#endif //PACE2024EXACT_ONE_VERTEX_PARTITION_MANAGER_H
