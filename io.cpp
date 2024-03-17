#include "io.h"
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <unordered_set>

Graph* readGraph(std::string graph_file) {
    //std::cout << "Reading graph..." << std::endl;
    std::string tmp; //for p and ocr in input format
    int n0; //neighbours = A, fixed partition
    int n1; //movable_nodes = B, free partition
    int m;

    std::ifstream file(graph_file);

    if (!file.is_open()) {
        //std::cout << "error: file not open" << std::endl;
    }

    std::string line;

    std::getline(file, line);

    while (line[0] == 'c') {
        std::getline(file, line);
    }

    std::stringstream ss(line);
    ss >> tmp; //p
    ss >> tmp; //ocr
    ss >> n0;
    ss >> n1;
    ss >> m;

    //std::cout << "n0: " << n0 << std::endl;
    //std::cout << "n1: " << n1 << std::endl;
    //std::cout << "m: " << m << std::endl;

    //initialize graph
    Graph* g = new Graph(n0, n1, m);
    //read adjacencies of the nodes in the graph file
    while (std::getline(file, line)) {
        if (line[0] == 'c' || line.empty()) {
            continue;
        }
        std::stringstream ss(line);

        int x;
        int y;

        ss >> x;
        ss >> y;
        g->addEdge(x-1, y-1);
    }
    file.close();

    //g->printGraph();
    g->sortNeighbours();
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
        outputFile << node->id + 1 << std::endl;
    }

    outputFile.close();
}

void writeGraphToBipartiteGraph(int n0, int n1, int m, std::vector<std::pair<int, int>> edges, std::string output){
    std::ofstream outputFile(output);

    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        return;
    }

    // Write header comments
    outputFile << "c first number: size neighbours (A)\n";
    outputFile << "c second number: size movable_nodes (B)\n";
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
    
    writeGraphToBipartiteGraph(n0, n1, edges_in_bipartite_graph.size(), edges_in_bipartite_graph, output); 

    file.close();
    return;
}

void readWeightedHyperGraph(std::string hypergraph_file, std::string output){
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
    int n0, n1, fmp;

    ss >> n0;
    ss >> n1;
    ss >> fmp;

    //std::cout << "n0: " << n0 << " n1: "<< n1 << std::endl;

    std::vector<std::pair<int, int>> edges_in_bipartite_graph;
    //read adjacencies of the nodes in the graph file
    int edge_index = 1;
    while (std::getline(file, line)) {
        // Skip comments
        if (line[0] == '%') {
            continue;
        }
    
        std::stringstream ss(line);
        int weight;
        ss >> weight;

        int edge_vertex;
        while (ss >> edge_vertex){
            //Need to offset edge_vertex indices by ammount of vertices in A partition (= n0)
            edges_in_bipartite_graph.push_back(std::make_pair(edge_index, n0 + edge_vertex));
        }          

        //Prep next edges
        edge_index++;
    }
    
    writeGraphToBipartiteGraph(n0, n1, edges_in_bipartite_graph.size(), edges_in_bipartite_graph, output); 

    file.close();
    return;
}

std::vector<std::pair<int, int>> generateBipartiteGraph(int sizeX, int sizeY) {
    std::vector<std::pair<int, int>> edges;

    // Create vectors to store nodes in each partition
    std::vector<int> partitionA(sizeX);
    std::vector<int> partitionB(sizeY);

    // Fill partition vectors with node indices
    std::iota(partitionA.begin(), partitionA.end(), 1);
    std::iota(partitionB.begin(), partitionB.end(), sizeX + 1);

    // Shuffle the partition vectors randomly
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(partitionA.begin(), partitionA.end(), gen);
    std::shuffle(partitionB.begin(), partitionB.end(), gen);

    // Create a uniform distribution for the number of edges per node
    std::uniform_int_distribution<> disEdges(1, std::min(sizeY, 5));

    // Generate random edges between partitions
    for (int x : partitionA) {
        // Generate a random number of edges for node x
        int numEdges = disEdges(gen);
        
        // Create a set to keep track of already chosen nodes in partition B
        std::unordered_set<int> chosenNodes;
        
        // Generate random edges for node x
        for (int i = 0; i < numEdges; ++i) {
            int y = partitionB[gen() % sizeY]; // Random node in partition B
            
            // Ensure uniqueness of edges
            while (chosenNodes.count(y) > 0) {
                y = partitionB[gen() % sizeY];
            }
            chosenNodes.insert(y);
            
            edges.push_back({x, y});
        }
    }

    std::sort(edges.begin(), edges.end());

    return edges;
}