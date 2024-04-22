#ifndef PERFORMANCE_H
#define PERFORMANCE_H

#include "graph.h"

void writeCSV(const std::string& filename, const std::vector<std::vector<std::string>>& data);
void calculatePerformance(std::string folderPath, std::string outputFile, int mode = 0, int method1 = 0, int method2 = 0, bool fast = true, bool almost = true);




#endif //PERFORMANCE_H