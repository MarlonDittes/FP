#ifndef PACE2024EXACT_PERFECT_PATTERN_H
#define PACE2024EXACT_PERFECT_PATTERN_H

#include <algorithm>
#include <cstring>

#include "definitions.h"
#include "macros.h"
#include "misc.h"

namespace CrossGuard {

    class PerfectPattern {
    private:
        u32 m_n = 0;
        AlignedVector<u8> m_matches;

    public:
        explicit PerfectPattern(u32 n) {
            m_n = n;
            m_matches.resize(n * n, 0);
        }

        inline void add(u32 i, u32 j) {
            ASSERT((i * m_n) + j < m_n * m_n);
            ASSERT(i < m_n);
            ASSERT(j < m_n);
            m_matches[i * m_n + j] = 1;
        }

        inline bool matches(u32 i, u32 j) {
            ASSERT((i * m_n) + j < m_n * m_n);
            ASSERT(i < m_n);
            ASSERT(j < m_n);
            return m_matches[i * m_n + j];
        }
    };

}

#endif //PACE2024EXACT_PERFECT_PATTERN_H
