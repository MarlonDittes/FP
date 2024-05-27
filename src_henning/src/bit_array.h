#ifndef PACE2024EXACT_BIT_ARRAY_H
#define PACE2024EXACT_BIT_ARRAY_H

#include <iostream>
#include <vector>

#include "definitions.h"
#include "macros.h"
#include "misc.h"

namespace CrossGuard {

    class BitArray {
    private:
        AlignedVector<u32> m_vec;

    public:
        BitArray() = default;

        void resize(u32 n){
            u32 vec_size = n / 32 + (n % 32 != 0);
            m_vec.resize(vec_size);
            std::fill(m_vec.begin(), m_vec.end(), 0);
        }

        void set(u32 i) {
            u32 idx = i / 32;
            u32 offset = i % 32;
            m_vec[idx] |= (1 << offset);
        }

        void clear(u32 i) {
            u32 idx = i / 32;
            u32 offset = i % 32;
            m_vec[idx] &= ~(1 << offset);
        }

        bool get(u32 i) const {
            u32 idx = i / 32;
            u32 offset = i % 32;
            return (m_vec[idx] >> offset) & 1u;
        }

    };
}

#endif //PACE2024EXACT_BIT_ARRAY_H
