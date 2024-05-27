#ifndef PACE2024EXACT_CANDIDATEMANAGER_H
#define PACE2024EXACT_CANDIDATEMANAGER_H

#include <vector>

#include "definitions.h"
#include "macros.h"
#include "misc.h"

namespace CrossGuard {

    /**
     * One candidate.
     */
    struct Candidate {
        u32 c; // the element
        u32 gain; // gain in the number of cuts
    };

    /**
     * Class to manage the candidates.
     */
    class CandidateManager {
    public:
        AlignedVector<Candidate> candidates;

        u32 max_median = 0;
        AlignedVector<std::pair<u32, u32>> median_candidate;

    public:
        /**
         * Default constructor.
         *
         * @param n Number of elements.
         */
        explicit CandidateManager(u32 n) {
            candidates.resize(n);
        }

        inline void initialize_median_vec(u32 n){
            median_candidate.resize(n);
            max_median = 0;
        }

        /**
         * Sorts the CandidateManager ascending.
         */
        inline void sort(u32 size) {
            AlignedVector<std::pair<u32, u32>> median_candidate_helper(median_candidate.size());
            AlignedVector<Candidate> candidates_helper(candidates.size());
            AlignedVector<u32> counting;

            // initialize counting
            counting.resize(max_median + 1);
            std::fill(counting.begin(), counting.end(), 0);

            // count
            for (u32 i = 0; i < size; ++i) {
                counting[median_candidate[i].second]++;
            }

            // prefix sum
            for (u32 i = 1; i <= max_median; ++i) {
                counting[i] += counting[i - 1];
            }

            // generate sorted output
            for (int i = (int) size - 1; i >= 0; i--) {
                median_candidate_helper[counting[median_candidate[i].second] - 1] = median_candidate[i];
                counting[median_candidate[i].second]--;
            }

            // generate solution
            for(u32 i = 0; i < median_candidate.size(); ++i){
                candidates_helper[i] = candidates[median_candidate_helper[i].first];
            }

            candidates.swap(candidates_helper);
            median_candidate.clear();
            median_candidate.shrink_to_fit();
        }

        inline void print(u32 size) {
            for (u32 i = 0; i < size; ++i) {
                std::cout << "(" << candidates[i].c << ", " << candidates[i].gain << "), ";
            }
            std::cout << std::endl;
        }

    };

}

#endif //PACE2024EXACT_CANDIDATEMANAGER_H
