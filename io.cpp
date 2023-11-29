#include "io.h"
#include <vector>

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

void writeHyperGraphToBipartiteGraph(int n0, int n1, int m, std::vector<std::pair<int, int>> edges, std::string output){
    std::ofstream outputFile(output);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return;
    }

    // Write header comments
    outputFile << "c first number: size X (A)\n";
    outputFile << "c second number: size Y (B)\n";
    outputFile << "c third number: number of edges (m)\n";
    outputFile << "p ocr " << n0 << " " << n1 << " " << m << "\n";

    // Write hypergraph edges
    for (const auto& edge : edges) {
        outputFile << edge.first << " " << edge.second << "\n";
    }

    outputFile.close();
    return;
}

void readHyperGraph(std::string hypergraph_file, std::string output){
    std::cout << "read hypergraph..." << std::endl;

    std::ifstream file(hypergraph_file);

    if (!file.is_open()) {
        std::cout << "error: file not open" << std::endl;
    }

    std::string line;
    std::getline(file, line);

    // Skipping comments
    while (line[0] == '%') {
        std::getline(file, line);
    }

    // First non comment line
    std::stringstream ss(line);
    int number_hyperedges, number_vertices;

    ss >> number_hyperedges;
    ss >> number_vertices;
    
    int n0 = number_hyperedges;
    int n1 = number_vertices;

    std::cout << "n0: " << n0 << " n1: "<< n1 << std::endl;

    std::vector<std::pair<int, int>> edges_in_bipartite_graph;
    //read adjacencies of the nodes in the graph file
    int edge_index = 1;
    while (std::getline(file, line)) {
        // Skip comments
        if (line[0] == '%') {
            continue;
        }
    
        std::stringstream ss(line);

        int edge_vertex;
        while (ss >> edge_vertex){
            //Need to offset edge_vertex indices by ammount of vertices in A partition (= n0)
            edges_in_bipartite_graph.push_back(std::make_pair(edge_index, n0 + edge_vertex));
        }          

        //Prep next edges
        edge_index++;
    }
    
    writeHyperGraphToBipartiteGraph(n0, n1, edges_in_bipartite_graph.size(), edges_in_bipartite_graph, output); 

    file.close();
    return;
}