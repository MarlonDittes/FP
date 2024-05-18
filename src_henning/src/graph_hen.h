#ifndef PACE2024EXACT_GRAPH_H
#define PACE2024EXACT_GRAPH_H

#include <cstdint>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <sstream>

#include "definitions.h"
#include "macros.h"
#include "misc.h"
#include "edge.h"

namespace CrossGuard {
    /**
     * Class to store the bipartite graph.
     */
    class Graph {
    public:
        u32 n_A = 0;
        u32 n_B = 0;
        u32 n_edges = 0;

        AlignedVector<AlignedVector<Edge>> adj_list;
        AlignedVector<AlignedVector<Edge>> adj_list_from_A;

        AlignedVector<u32> medians;
        AlignedVector<u64> adj_hash;

        bool is_finalized = false;


        /**
         * Reads a graph from a file.
         *
         * @param file_path Path to the file.
         */
        explicit Graph(const std::string &file_path) {
            std::ios_base::sync_with_stdio(false);
            std::ifstream file(file_path);

            std::string line(64, ' ');
            if (file.is_open()) {
                while (std::getline(file, line)) {
                    if (line[0] == 'c') { continue; }
                    if (line.back() == '\n' || line.back() == '\r') { line.pop_back(); }

                    assert(line[0] == 'p');
                    assert(line[1] == ' ');
                    assert(line[2] == 'o');
                    assert(line[3] == 'c');
                    assert(line[4] == 'r');
                    assert(line[5] == ' ');

                    n_A = 0;
                    n_B = 0;
                    n_edges = 0;
                    size_t i = 6;

                    // read in n_A
                    while (line[i] != ' ') {
                        n_A = n_A * 10 + line[i] - '0';
                        i += 1;
                    }
                    i += 1;

                    // read in n_B
                    while (line[i] != ' ') {
                        n_B = n_B * 10 + line[i] - '0';
                        i += 1;
                    }

                    // adjacency lists
                    adj_list_from_A.resize(n_A);
                    adj_list.resize(n_B);
                    break;
                }

                while (std::getline(file, line)) {
                    if (line[0] == 'c') { continue; }
                    if (line.back() == '\n' || line.back() == '\r') { line.pop_back(); }

                    u32 a = 0, b = 0;
                    u64 i = 0;

                    // read in vertex a
                    while (line[i] != ' ') {
                        a = a * 10 + line[i] - '0';
                        i += 1;
                    }
                    a -= 1;
                    i += 1;

                    // read in vertex b
                    while (i < line.size()) {
                        b = b * 10 + line[i] - '0';
                        i += 1;
                    }

                    b -= (n_A + 1);
                    add_edge(a, b, 1);
                }
                file.close();

                // sort each neighborhood
                finalize();
            } else {
                std::cout << "Could not open file " << file_path << " !" << std::endl;
            }
        }

        /**
         * Initializes the graph.
         *
         * @param n_A Number of fixed vertices.
         * @param n_B Number of movable vertices.
         */
        Graph(u32 n_A, u32 n_B) {
            this->n_A = n_A;
            this->n_B = n_B;
            n_edges = 0;

            adj_list.resize(n_B);
            adj_list_from_A.resize(n_A);
        }

        Graph() {
            n_A = 0;
            n_B = 0;
            n_edges = 0;
        }

        /**
         * Adds an edge to the graph. Note that a must be in the range
         * {0, ..., n_A - 1} and b must be in the range {0, ..., n_B - 1}.
         *
         * @param a Vertex a (fixed vertex).
         * @param b Vertex b (movable vertex).
         */
        inline void add_edge(u32 a, u32 b, u32 weight) {
            ASSERT(a < n_A);
            ASSERT(b < n_B);
            ASSERT(a < (u32) adj_list_from_A.size());
            ASSERT(b < (u32) adj_list.size());

            n_edges += 1;
            adj_list[b].push_back({a, weight});
            adj_list_from_A[a].push_back({b, weight});

            ASSERT(no_duplicates(adj_list[b]));
            ASSERT(no_duplicates(adj_list_from_A[a]));
        }

        /**
         * Finalized the graph.
         */
        inline void finalize() {
            // sort each neighborhood
            sort_neighborhood();

            // determine the median of each B vertex
            calc_medians();

            // determine the hash for each neighborhood
            calc_adj_hashes();

            is_finalized = true;
        }

        /**
         * Determines the number of cuts. Note that the vector should already be
         * shifted in the range [0, ..., m_n_B - 1].
         *
         * @param permutation The m_permutation.
         * @return Number of cuts.
         */
        inline u32 determine_n_cuts(const AlignedVector<u32> &permutation) {
            ASSERT((u32) permutation.size() == n_B);
            ASSERT(no_duplicates(permutation));

            u32 n_cuts = 0;
            for (u32 i = 0; i < n_B; ++i) {
                for (u32 j = i + 1; j < n_B; ++j) {
                    u32 b1 = permutation[i];
                    u32 b2 = permutation[j];

                    // loop through the edges
                    for (size_t k = 0; k < adj_list[b1].size(); ++k) {
                        for (size_t l = 0; l < adj_list[b2].size(); ++l) {
                            Edge a1 = adj_list[b1][k];
                            Edge a2 = adj_list[b2][l];

                            n_cuts += (a2.vertex < a1.vertex) * (a1.weight * a2.weight);
                        }
                    }
                }
            }

            return n_cuts;
        }

        /**
         * Prints the graph as an adjacency list. Only prints the adjacency of
         * movable nodes.
         */
        inline void print() const {
            for (u32 i = 0; i < n_B; ++i) {
                std::cout << i << ": ";
                CrossGuard::print(adj_list[i]);
            }
        }

        /**
         * Returns the median solution.
         *
         * @return The median solution.
         */
        inline AlignedVector<unsigned int> get_median_solution(){
            ASSERT(is_finalized);

            AlignedVector<std::pair<u32, u32>> temp(n_B);
            for(u32 i = 0; i < n_B; ++i){
                u32 vertex_b = i;
                u32 median = medians[i];
                temp[i] = {vertex_b, median};
            }

            std::sort(temp.begin(), temp.end(), [](auto &left, auto &right) {
                return left.second < right.second;
            });

            AlignedVector<u32> res(n_B);
            for(u32 i = 0; i < n_B; ++i){
                res[i] = temp[i].first;
            }

            return res;
        }

    private:
        /**
         * Sorts each neighborhood ascending.
         */
        inline void sort_neighborhood() {
            // sort each neighborhood
            for (u32 i = 0; i < n_B; ++i) {
                std::sort(adj_list[i].begin(), adj_list[i].end());
            }

            for (u32 i = 0; i < n_A; ++i) {
                std::sort(adj_list_from_A[i].begin(), adj_list_from_A[i].end());
            }

            for (u32 i = 0; i < n_B; ++i) {
                ASSERT(std::is_sorted(adj_list[i].begin(), adj_list[i].end()));
            }

            for (u32 i = 0; i < n_A; ++i) {
                ASSERT(std::is_sorted(adj_list_from_A[i].begin(), adj_list_from_A[i].end()));
            }
        }

        /**
         * Calculates the median of each vertex in B.
         */
        inline void calc_medians(){
            medians.resize(n_B);

            for(u32 vertex_b = 0; vertex_b < n_B; ++vertex_b){
                u32 median = 0;
                if(!adj_list[vertex_b].empty()){
                    u32 size = adj_list[vertex_b].size();
                    median = adj_list[vertex_b][size / 2].vertex;
                }
                medians[vertex_b] = median;
            }
        }

        /**
         * Calculates the hashes of each neighborhood.
         */
        inline void calc_adj_hashes(){
            adj_hash.resize(n_B);

            for(u32 i = 0; i < n_B; ++i){
                adj_hash[i] = hash(adj_list[i]);
            }
        }

    };


}

#endif //PACE2024EXACT_GRAPH_H
