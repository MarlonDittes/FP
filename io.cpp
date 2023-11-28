#include "io.h"

Graph* readGraph(std::string graph_file) {

    std::cout << "read graph..." << std::endl;
    std::string tmp; //for p and ocr in input format
    int size_X; //X = A, fixed partition
    int size_Y; //Y = B, free partition
    int edge_number;

    std::ifstream file(graph_file);

    if (!file.is_open()) {
        std::cout << "error: file not open" << std::endl;
    }

    std::string line;

    std::getline(file, line);

    while (line[0] == 'c') {
        std::cout << "skip" << std::endl;
        std::getline(file, line);
    }

    std::stringstream ss(line);
    ss >> tmp; //p
    ss >> tmp; //ocr
    std::cout << "tmp: " << tmp << std::endl;
    ss >> size_X;
    ss >> size_Y;
    ss >> edge_number;

    std::cout << "size_X: " << size_X << std::endl;
    std::cout << "size_Y: " << size_Y << std::endl;
    std::cout << "edge_number: " << edge_number << std::endl;

    //initialize graph
    Graph* g = new Graph(size_X, size_Y, edge_number);
    std::cout << "test" << std::endl;
    //read adjacencies of the nodes in the graph file
    while (std::getline(file, line)) {
        std::cout << "while..." << std::endl;

        if (line[0] == 'c') {
            continue;
        }
        std::stringstream ss(line);

        int x;
        int y;

        ss >> x;
        ss >> y;
        std::cout << "addEdge: " << x << " - " << y << std::endl;
        g->addEdge(x, y);
    }
    std::cout << "close" << std::endl;
    file.close();

    //g.printGraph();
    return g;
}

void outputOrder(std::vector<Node*> order, std::string output) {
    std::ofstream outputFile(output);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return;
    }

    for (const auto& node : order) {
        // Print the node id to the file
        outputFile << node->id << std::endl;
    }

    outputFile.close();
}