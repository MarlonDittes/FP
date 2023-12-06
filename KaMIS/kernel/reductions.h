/******************************************************************************
* reductions.h
*
*****************************************************************************/

#ifndef REDUCTIONS_H
#define REDUCTIONS_H

// local includes
#include "definitions.h"
#include "fast_set.h"
#include "data_structure/sized_vector.h"
#include "data_structure/dynamic_graph.h"

// system includes
#include <vector>
#include <memory>
#include <array>

class branch_and_reduce_algorithm;

enum reduction_type {twin_reduction};
constexpr size_t REDUCTION_NUM = 1;

class vertex_marker {
private:
	sized_vector<NodeID> current;
	sized_vector<NodeID> next;
	fast_set added_vertices;

public:
	vertex_marker(size_t size) : current(size), next(size), added_vertices(size) {};

	void add(NodeID vertex) {
		if (!added_vertices.get(vertex)) {
			next.push_back(vertex);
			added_vertices.add(vertex);
		}
	}

	void get_next() {
		if (next.size() != 0) {
			current.swap(next);
			clear_next();
		}
	}

	void clear_next() {
		next.clear();
		added_vertices.clear();
	}

	void fill_current_ascending(size_t n) {
		current.clear();
		for (size_t i = 0; i < n; i++) {
			current.push_back(i);
		}
	}

	NodeID current_vertex(size_t index) {
		return current[index];
	}

	size_t current_size() {
		return current.size();
	}
};

struct general_reduction {
	general_reduction(size_t n) : marker(n) {}
	virtual ~general_reduction() {}
	virtual general_reduction* clone() const = 0;

	virtual reduction_type get_reduction_type() const = 0;
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) = 0;
	virtual void restore(branch_and_reduce_algorithm* br_alg) {}
	virtual void apply(branch_and_reduce_algorithm* br_alg) {}

	bool has_run = false;
	vertex_marker marker;
};

struct twin_reduction : public general_reduction {
	twin_reduction(size_t n) : general_reduction(n) {}
	~twin_reduction() {}
	virtual twin_reduction* clone() const final { return new twin_reduction(*this); }

	virtual reduction_type get_reduction_type() const final { return reduction_type::twin; }
	virtual bool reduce(branch_and_reduce_algorithm* br_alg) final;
	virtual void restore(branch_and_reduce_algorithm* br_alg) final;
	virtual void apply(branch_and_reduce_algorithm* br_alg) final;

private:
	struct restore_data {
		NodeID main;
		NodeID twin;
	};

	void fold(branch_and_reduce_algorithm* br_alg, NodeID main, NodeID twin);

	std::vector<restore_data> restore_vec;
};

struct reduction_ptr {
	general_reduction* reduction = nullptr;

	reduction_ptr() = default;

	~reduction_ptr() {
		release();
	}

	reduction_ptr(general_reduction* reduction) : reduction(reduction) {};

	reduction_ptr(const reduction_ptr& other) : reduction(other.reduction->clone()) {};

	reduction_ptr& operator=(const reduction_ptr& other) {
		release();
		reduction = other.reduction->clone();
		return *this;
	};

	reduction_ptr(reduction_ptr&& other) : reduction(std::move(other.reduction)) {
		other.reduction = nullptr;
	};

	reduction_ptr& operator=(reduction_ptr&& other) {
		reduction = std::move(other.reduction);
		other.reduction = nullptr;
		return *this;
	};

	general_reduction* operator->() const {
		return reduction;
	}

	void release() {
		if (reduction) {
			delete reduction;
			reduction = nullptr;
		}
	};
};

template<class Last>
void make_reduction_vector_helper(std::vector<reduction_ptr>& vec, size_t n) {
	vec.emplace_back(new Last(n));
};

template<class First, class Second, class ...Redus>
void make_reduction_vector_helper(std::vector<reduction_ptr>& vec, size_t n) {
	vec.emplace_back(new First(n));
	make_reduction_vector_helper<Second, Redus...>(vec, n);
};

template<class ...Redus>
std::vector<reduction_ptr> make_reduction_vector(size_t n) {
	std::vector<reduction_ptr> vec;
	make_reduction_vector_helper<Redus...>(vec, n);
	return vec;
};


/*template<class ...Redus>
std::vector<std::unique_ptr<general_reduction>> make_reduction_vector(size_t n) {
std::vector<std::unique_ptr<general_reduction>> vec;
(vec.push_back(std::make_unique<Redus>(n)), ...);
return vec;
}

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class Last>
void make_reduction_vector_helper(std::vector<std::unique_ptr<general_reduction>>& vec, size_t n) {
vec.push_back(make_unique<Last>(n));
};

template<class First, class Second, class ...Redus>
void make_reduction_vector_helper(std::vector<std::unique_ptr<general_reduction>>& vec, size_t n) {
vec.push_back(make_unique<First>(n));
make_reduction_vector_helper<Second, Redus...>(vec, n);
};

template<class ...Redus>
std::vector<std::unique_ptr<general_reduction>> make_reduction_vector(size_t n) {
std::vector<std::unique_ptr<general_reduction>> vec;
make_reduction_vector_helper<Redus...>(vec, n);
return vec;
};*/


#endif //REDUCTIONS_H
