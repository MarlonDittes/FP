#include "stack.h"
#include <iostream>
#include <stdexcept>

// Push element onto stack
void Stack::push(int data) {
    data_vec.push_back(data);
}

// Check if the stack is empty
bool Stack::isEmpty() const {
    return data_vec.size() == 0;
}

// Pop element from stack
int Stack::pop() {
    if (isEmpty()) {
        throw std::out_of_range("Stack Underflow");
    }

    int to_return_value = data_vec[data_vec.size() - 1];
    data_vec.pop_back();
    return to_return_value;
}

// Peek at the top element of stack
int Stack::peek() const {
    if (isEmpty()) {
        throw std::out_of_range("Stack is Empty");
    }

    return data_vec[data_vec.size() - 1];
}

