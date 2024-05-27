#ifndef PACE2024EXACT_MISC_H
#define PACE2024EXACT_MISC_H

#include <vector>
#include <iostream>
#include <cstdint>
#include <fstream>
#include <chrono>
#include <sys/stat.h>
#include <unistd.h>
#include <string>

#include "definitions.h"
#include "macros.h"
#include "edge.h"

namespace CrossGuard {

/**
 * Checks if a file exists.
 *
 * @param file_path Path to the file.
 * @return True if file exists, False else.
 */
    bool file_exists(const std::string &file_path);

    /**
 * Prints the content of a vector in a Numpy-style way. The '<<' operator must
 * be defined for T.
 *
 * @tparam T Type of the vector.
 * @param vec The vector.
 */
    template<typename T>
    void print(const AlignedVector<T> &vec) {
        if (vec.empty()) {
            std::cout << "[]" << std::endl;
            return;
        }

        std::cout << "[";
        for (size_t i = 0; i < vec.size() - 1; ++i) {
            std::cout << vec[i] << ", ";
        }
        std::cout << vec.back() << "]" << std::endl;
    }

/**
 * Returns the content of a vector in a Numpy-style way. The '<<' operator must
 * be defined for T.
 *
 * @tparam T Type of the vector.
 * @param vec The vector.
 */
    template<typename T>
    std::string to_string(const AlignedVector<T> &vec) {
        if (vec.empty()) {
            return "[]";
        }
        std::string s = "[";

        for (size_t i = 0; i < vec.size() - 1; ++i) {
            s += std::to_string(vec[i]) + ", ";
        }
        return s + std::to_string(vec.back()) + "]\n";
    }

/**
 * Prints the content of a vector in a Numpy-style way. The '<<' operator must
 * be defined for T.
 *
 * @tparam T Type of the vector.
 * @param vec The vector.
 * @param size The size of the vector.
 */
    template<typename T>
    void print(const AlignedVector<T> &vec, size_t size) {
        if (size == 0) {
            std::cout << "[]" << std::endl;
            return;
        }

        std::cout << "[";
        for (size_t i = 0; i < size - 1; ++i) {
            std::cout << vec[i] << ", ";
        }
        std::cout << vec[size - 1] << "]" << std::endl;
    }

/**
 * Determines if a vector contains a duplicate item. The '==' operator must be
 * defined for T.
 *
 * @tparam T Type of the vector.
 * @param vec The vector.
 * @return True if no duplicates are present, false else.
 */
    template<typename T>
    bool no_duplicates(const AlignedVector<T> &vec) {
        for (size_t i = 0; i < vec.size(); ++i) {
            for (size_t j = i + 1; j < vec.size(); ++j) {
                if (vec[i] == vec[j]) {
                    return false;
                }
            }
        }
        return true;
    }

/**
 * Specialized function that determines if each vector inside the vector is
 * unique.
 *
 * @param vec Vector of vectors.
 * @return True if no duplicates are present, false else.
 */
    bool no_duplicates(const AlignedVector<AlignedVector<int>> &vec);

/**
 * Computes the factorial of a number.
 *
 * @tparam T Type of the number.
 * @param n The number.
 * @return n!
 */
    template<typename T>
    T fac(T n) {
        int prod = 1;
        for (T i = 1; i <= n; ++i) {
            prod *= i;
        }
        return prod;
    }

/**
 * Reads a m_solution from a .sol file.
 *
 * @param file_path Path to the file.
 * @param shift Shifts the read numbers to the left (x -= shift).
 * @return Vector containing the m_solution.
 */
    AlignedVector<unsigned int> read_solution(const std::string &file_path, unsigned int shift);

/**
 * Write a m_solution vector to a file.
 *
 * @param solution Vector holding the m_solution.
 * @param file_path Path to the file.
 */
    void write_solution(const AlignedVector<unsigned int> &solution, const std::string &file_path);

    /**
     * Returns the difference between two time points in seconds.
     *
     * @param sp Start time point.
     * @param ep End time point.
     * @return Difference in seconds.
     */
    double get_seconds(std::chrono::steady_clock::time_point sp,
                       std::chrono::steady_clock::time_point ep);

    template<typename T>
    bool exists(const AlignedVector<T> &vec, const T &e) {
        for (auto &x: vec) {
            if (x == e) {
                return true;
            }
        }
        return false;
    }

    template<typename T>
    bool equal(const AlignedVector<T> &vec1, const AlignedVector<T> &vec2) {
        if (vec1.size() != vec2.size()) {
            return false;
        }

        for (size_t i = 0; i < vec1.size(); ++i) {
            if (vec1[i] != vec2[i]) {
                return false;
            }
        }
        return true;
    }

    template<typename T>
    T min(const AlignedVector<T> &vec) {
        T m = vec[0];

        for (auto &x: vec) {
            m = std::min(m, x);
        }

        return m;
    }

    template<typename T>
    T max(const AlignedVector<T> &vec) {
        T m = vec[0];

        for (auto &x: vec) {
            m = std::max(m, x);
        }

        return m;
    }

    template<typename T>
    u64 hash(const AlignedVector<T> &vec) {
        std::size_t seed = vec.size();
        for(auto x : vec) {
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    u64 hash(const AlignedVector<Edge> &vec);

}

#endif //PACE2024EXACT_MISC_H
