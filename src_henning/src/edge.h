#ifndef PACE2024EXACT_EDGE_H
#define PACE2024EXACT_EDGE_H

#include <iostream>

#include "definitions.h"

namespace CrossGuard {
    /**
     * Represents one edge.
     */
    class Edge {
    public:
        u32 vertex;
        u32 weight;
    };

    bool operator==(const Edge &lhs, const Edge &rhs);

    bool operator<=(const Edge &lhs, const Edge &rhs);

    bool operator<(const Edge &lhs, const Edge &rhs);

    std::ostream &operator<<(std::ostream &os, const Edge &obj);

}

#endif //PACE2024EXACT_EDGE_H
