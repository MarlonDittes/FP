#include <iostream>
#include <filesystem>
#include <fstream>
#include <chrono>

#include "graph.h"
#include "reductions.h"
#include "io.h"

namespace fs = std::filesystem;

void writeCSV(const std::string& filename, const std::vector<std::vector<std::string>>& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Iterate through the data and write each row to the file
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            // Add a comma if it's not the last column
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n"; // Add a newline after each row
    }

    file.close();
}

void calculatePerformance(std::string folderPath, std::string outputFile, int mode, int method1, int method2, bool fast, bool almost){
    /* Modes:   1 BruteForce 
                2 Barycenter
                3 Median
                4 BranchAndReduce
    */

    try {
        std::vector<std::vector<std::string>> data = {{"graph", "duration[μs]", "crossings",
                                                         "ZeroEdge_usages", "Complete_usages", "ZeroCrossings_usages",
                                                         "Twin_usages", "AlmostTwin_usages"}};
        // Iterating through all graphs in folderPath
        int count = 0;
        for (const auto& entry : fs::directory_iterator(folderPath)) {
            std::cout << count << std::endl;
            //if (count == 50){
            //    break;
            //}
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

                //Define early so switch works, used to measure elapsed time
                auto start = std::chrono::high_resolution_clock::now();
                std::pair<std::vector<Node *>, long> result;
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

                // Setup which reductions should be used
                std::vector<general_reduction*> reductions;
                reductions.push_back(new ZeroEdge_reduction);
                reductions.push_back(new Complete_reduction);
                reductions.push_back(new ZeroCrossings_reduction);
                reductions.push_back(new Twins_reduction);
                if (almost){
                    reductions.push_back(new AlmostTwin_reduction);
                }
                //reductions.push_back(new Domination_reduction); 

                switch (mode) {
                    // 1 BruteForce 
                    case 1:
                        start = std::chrono::high_resolution_clock::now();
                        result = bruteForce(g);
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
                        break;
                    // 2 Barycenter
                    case 2:
                        start = std::chrono::high_resolution_clock::now();
                        g->BarycenterHeuristicMarlon();
                        result = std::make_pair(g->getOrderNodes(), g->countCrossings());
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
                        break;
                    // 3 Median
                    case 3:
                        start = std::chrono::high_resolution_clock::now();
                        g->MedianHeuristic();
                        result = std::make_pair(g->getOrderNodes(), g->countCrossings());
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
                        break;
                    // 4 BranchAndReduce
                    case 4:
                        start = std::chrono::high_resolution_clock::now();
                        result = BranchAndReduce(g, reductions, method1, method2, fast);
                        stop = std::chrono::high_resolution_clock::now();
                        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
                        if (fast){
                            result.second = g->countCrossingsMarlon();
                        }
                        break;

                    default:
                        std::cout << "Select one of the following modes:" << std::endl;
                        std::cout << "1 BruteForce" << std::endl;
                        std::cout << "2 Barycenter" << std::endl;
                        std::cout << "3 Median" << std::endl;
                        std::cout << "4 BranchAndReduce" << std::endl;
                        return;
                }
                std::vector<std::string> newData = {testname, std::to_string(duration.count()), std::to_string(result.second),
                                                    std::to_string(reductions[0]->usage_count), std::to_string(reductions[1]->usage_count), std::to_string(reductions[2]->usage_count),
                                                    std::to_string(reductions[3]->usage_count), std::to_string(reductions[4]->usage_count)};
                data.push_back(newData);
            }
            count++;
        }
        writeCSV(outputFile, data);
    } catch (const fs::filesystem_error& ex) {
        std::cerr << ex.what() << std::endl;
    }

    
    
}