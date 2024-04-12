#ifndef REDUCTIONS_H
#define REDUCTIONS_H
#include <vector>
#include "graph.h"

class Graph;

enum reduction_type { ZeroEdge, Complete, ZeroCrossings, Twins, AlmostTwins, Domination };

struct general_reduction {
    general_reduction() {}
    virtual reduction_type get_reduction_type() const = 0;
	virtual int reduce(Graph* g) = 0;
    virtual bool apply(Graph* g, int twins_count) {return false;};

    int usage_count = 0;
};

//naive reduction = graph is optimal and will be no more reduced, no apply function
struct ZeroEdge_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::ZeroEdge; }
    virtual int reduce(Graph* g) override;
};

//naive reduction
struct Complete_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::Complete; }
    virtual int reduce(Graph* g) override;
};

//naive reduction
struct ZeroCrossings_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::ZeroCrossings; }
    virtual int reduce(Graph* g) override;
}; 

struct Twins_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::Twins; }
    virtual int reduce(Graph* g) override;
    virtual bool apply(Graph* g, int twins_count) override;

    private:
	    struct restore_data {
		    int main; //save nodeID
		    int twin;
	    };

	std::vector<restore_data> restore_vec; //saves twins as pairs of main and twin in a vector
};

struct AlmostTwin_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::AlmostTwins; }
    virtual int reduce(Graph* g) override;
    virtual bool apply(Graph* g, int twins_count) override;

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
    virtual int reduce(Graph* g) override;
    virtual bool apply(Graph* g, int twins_count) override;

    private:
	    struct restore_data {
		    int main; //save nodeID ->maybe needs to be changed to Nodes
		    int part;
            bool side; //which side to place twin s. t. optimal: 0 = left, 1 = right
            //save additional information for apply function:
            int domination_size; //number of identical neighbours
            int start; //id of first identical neighbour
	    };

	std::vector<restore_data> restore_vec; //saves twins as pairs of main and twin in a vector
};

#endif //REDUCTIONS_H