#include "stack.h"
#include <iostream>
#include <stdexcept>

Stack::Stack() : top(nullptr) {}

    // Destructor to deallocate memory
Stack::~Stack() {
    while (top != nullptr) {
        Stack_Node* temp = top;
        top = top->next;
        delete temp;
    }
}

    // Push element onto stack
void Stack::push(int data) {
    Stack_Node* newNode = new Stack_Node;
    newNode->data = data;
    newNode->next = top;
    top = newNode;
}

// Check if the stack is empty
bool Stack::isEmpty() const {
    return top == nullptr;
}

// Pop element from stack
int Stack::pop() {
    if (isEmpty()) {
        throw std::out_of_range("Stack Underflow");
    }

    Stack_Node* temp = top;
    int poppedValue = top->data;
    top = top->next;
    delete temp;
    return poppedValue;
}

// Peek at the top element of stack
int Stack::peek() const {
    if (isEmpty()) {
        throw std::out_of_range("Stack is Empty");
    }

    return top->data;
}

