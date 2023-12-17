#ifndef REDUCTIONS_H
#define REDUCTIONS_H

#include <vector>

enum reduction_type { twin_reduction, cheap_reduction };

struct general_reduction {
    general_reduction() {}
    virtual reduction_type get_reduction_type() const = 0;
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) = 0;
};

struct zeroCrossings_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::cheap_reduction; }
    void reduce() {};
};


struct reductions {
    std::vector<general_reduction> reduction_stack;

    reductions() {
        reduction_stack = std::vector<general_reduction>();
    }
};


/*struct cheap_reduction : public general_reduction {
    virtual reduction_type get_reduction_type() const final { return reduction_type::cheap_reduction; }
    void reduce() {};
};*/

#endif //REDUCTIONS_H