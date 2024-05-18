#include "edge.h"

namespace CrossGuard {


    bool operator==(const CrossGuard::Edge &lhs, const CrossGuard::Edge &rhs) {
        return lhs.vertex == rhs.vertex;
    }

    bool operator<=(const CrossGuard::Edge &lhs, const CrossGuard::Edge &rhs) {
        return lhs.vertex <= rhs.vertex;
    }

    bool operator<(const CrossGuard::Edge &lhs, const CrossGuard::Edge &rhs) {
        return lhs.vertex < rhs.vertex;
    }

    std::ostream &operator<<(std::ostream &os, const CrossGuard::Edge &obj) {
        os << "(" << obj.vertex << ", " << obj.weight << ")";
        return os;
    }

}
