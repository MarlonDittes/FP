PACE CHALLENGE 2024 - Student Submission: heiCross

DESCRIPTION
The one-sided crossing minimization problem (OCM) asks for an order of one layer of a two-layer bipartite graph, 
where the other layer of nodes is fixed, such that the number of edge crossings is minimized. 
Our heiCross solver employs a branch and reduce approach specifically designed for the heuristic track of the challenge. 
It starts by partitioning the graph via articulation points. Then, three reduction rules and further partitioning recursivly 
reduce the size of the graph. If this approach can no longer be executed, we use branching to build a heuristic solution. 
When the partitions reach a specific size, an exact solver is used. 
For a detailed explanation of the heiCross solver, please refer to the 'PACE2024_heiCross.pdf' included in this repository.

DOI
https://doi.org/10.5281/zenodo.12507916

INSTALLATION
1. Clone the repository: 
    git clone https://github.com/MarlonDittes/heiCross.git
2. Navigate to the project directory: 
    cd heiCross
3. Install the Boost Library (required for the algorithm):
    On Ubuntu: sudo apt-get install libboost-all-dev
    On macOS: brew install boost
    On Windows: Download and install from the Boost Webiste https://www.boost.org/
4. Use the standard cmake build process:
    mkdir build
    cd build 
    cmake    
    make 
5. Run with:
    ./heiCross < InputFile > OutputFile
    If the OutputFile is omitted, the output is printed directly to the Console.

REQUIREMENTS
Boost Library https://www.boost.org/

AUTHORS
Marlon Dittes @MarlonDittes
Alvaro Garmendia @alvarogarmen
Tomer Haham @TomerHaham
Shai Peretz @Shaip161
Antonie Wagner @Tonyloni
Henning Woydt @HenningWoydt

CO-AUTHORS
Ernestine Großmann
Henrik Reinstädtler
Adil Chhabra
