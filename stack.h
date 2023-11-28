#pragma once
#include <iostream>
#include <stdexcept>

class Stack {
    struct Stack_Node {
        int data;
        Stack_Node* next;
    };

    Stack_Node* top;

public:
    
    Stack();

    // Destructor to deallocate memory
    ~Stack();

    // Push element onto stack
    void push(int data);

    // Check if the stack is empty
    bool isEmpty() const;

    // Pop element from stack
    int pop();

    // Peek at the top element of stack
    int peek() const;
};


