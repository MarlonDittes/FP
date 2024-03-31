#ifndef REDUCTIONS_H
#define REDUCTIONS_H
#include <vector>
#include "graph.h"

class Graph;

enum reduction_type { ZeroEdge, Complete, ZeroCrossings, Twins, AlmostTwins, Domination };
constexpr int BRUTE_CUTOFF = 10;

struct general_reduction {
    general_reduction() {}
    virtual reduction_type get_reduction_type() const = 0;
	virtual bool reduce(Graph* g) = 0;
    virtual void apply(Graph* g, int n = 0) {}
};

//naive reduction = graph is optimal and will be no more reduced, no apply function
struct ZeroEdge_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::ZeroEdge; }
    virtual bool reduce(Graph* g) override;
};

//naive reduction
struct Complete_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::Complete; }
    virtual bool reduce(Graph* g) override;
};

//naive reduction
struct ZeroCrossings_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::ZeroCrossings; }
    virtual bool reduce(Graph* g) override;
}; 

struct Twins_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::Twins; }
    virtual bool reduce(Graph* g) override;
    virtual void apply(Graph* g, int twinsCount) override;

    private:
	    struct restore_data {
		    int main; //save nodeID
		    int twin;
	    };

	std::vector<restore_data> restore_vec; //saves twins as pairs of main and twin in a vector
};

struct AlmostTwin_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::AlmostTwins; }
    virtual bool reduce(Graph* g) override;

    private:
	    struct restore_data {
		    int main; //save nodeID ->maybe needs to be changed to Nodes
		    int twin;
            bool side; //which side to place twin s. t. optimal: 0 = left, 1 = right
	    };

	std::vector<restore_data> restore_vec; //saves twins as pairs of main and twin in a vector
};

struct Domination_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::Domination; }
    virtual bool reduce(Graph* g) override;

    private:
	    struct restore_data {
		    int main; //save nodeID ->maybe needs to be changed to Nodes
		    int part;
            bool side; //which side to place twin s. t. optimal: 0 = left, 1 = right
	    };

	std::vector<restore_data> restore_vec; //saves twins as pairs of main and twin in a vector
};

#endif //REDUCTIONS_H