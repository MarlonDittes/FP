#include <iostream>
#include <filesystem>
#include <fstream>
#include <chrono>

#include "graph.h"
#include "io.h"

namespace fs = std::filesystem;

void calculatePerformance(std::string folderPath, std::string outputFile, int mode){
    /* Modes:   1 BruteForce 
                2 Greedy
                3 Baricenter
                4 Baricenter + Greedy
                5 BranchAndReduce
    */

    try {
        // Create file to save performance data in
        std::ofstream outFile(outputFile);
        if (!outFile.is_open()) {
            std::cerr << "Error opening "<< outputFile << std::endl;
            return;
        }
        // Iterating through all graphs in folderPath
        for (const auto& entry : fs::directory_iterator(folderPath)) {
            if (fs::is_regular_file(entry.path())) {
                size_t lastSlashPos = entry.path().string().find_last_of('/');
                std::string testname;
                if (lastSlashPos != std::string::npos) {
                    testname = entry.path().string().substr(lastSlashPos + 1);
                }

                //std::cout << "Testing on " << testname << std::endl;
                // Test performance on one graph
                Graph* g = readGraph(entry.path().string());
                g->sortNeighbours();

                //Define early so switch works
                auto start = std::chrono::high_resolution_clock::now();
                std::pair<std::vector<Node *>, long> result;
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

                std::vector<Reduction> reductions;

                switch (mode) {
                    // 1 BruteForce 
                    case 1:
                        start = std::chrono::high_resolution_clock::now();
                        result = bruteForce(g);
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

                        outFile << testname << " " << duration.count() << " " << result.second << std::endl;
                        break;
                    // 2 Greedy
                    case 2:
                        start = std::chrono::high_resolution_clock::now();
                        result = g->Greedy();
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

                        outFile << testname << " " << duration.count() << " " << result.second << std::endl;
                        break;
                    // 3 Baricenter
                    case 3:
                        start = std::chrono::high_resolution_clock::now();
                        g->Median_Heuristic();
                        result = std::make_pair(g->getOrderNodes(), g->countCrossings());
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

                        outFile << testname << " " << duration.count() << " " << result.second << std::endl;
                        break;
                    // 4 Baricenter + Greedy
                    case 4:
                        start = std::chrono::high_resolution_clock::now();
                        g->Median_Heuristic();
                        result = g->Greedy();
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                        
                        outFile << testname << " " << duration.count() << " " << result.second << std::endl;
                        break;
                    // 5 BranchAndReduce
                    case 5:
                        start = std::chrono::high_resolution_clock::now();
                        result = BranchAndReduce(g, reductions);
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
                        
                        outFile << testname << " " << duration.count() << " " << result.second << std::endl;
                        break;

                    default:
                        std::cout << "Select one of the following modes:" << std::endl;
                        std::cout << "1 BruteForce" << std::endl;
                        std::cout << "2 Greedy" << std::endl;
                        std::cout << "3 Baricenter" << std::endl;
                        std::cout << "4 Baricenter + Greedy" << std::endl;
                        return;
                }
            }
        }
        outFile.close();
    } catch (const fs::filesystem_error& ex) {
        std::cerr << ex.what() << std::endl;
    }

    
    
}