#ifndef PACE2024EXACT_DEFINITIONS_H
#define PACE2024EXACT_DEFINITIONS_H

#include <cstdint>
#include <cstdlib>
#include <vector>
#include <boost/align.hpp>

namespace CrossGuard {
    /**
     * All data is 64 byte aligned.
     */
    template<typename T>
    using AlignedVector = std::vector<T, boost::alignment::aligned_allocator<T, 64>>;
    // using AlignedVector = std::vector<T>;

    typedef int8_t s8;
    typedef int16_t s16;
    typedef int32_t s32;
    typedef int64_t s64;

    typedef u_int8_t u8;
    typedef u_int16_t u16;
    typedef u_int32_t u32;
    typedef u_int64_t u64;

    typedef float f32;
    typedef double f64;
}

#endif //PACE2024EXACT_DEFINITIONS_H
