PACE CHALLENGE 2024 - Team: HeiCross

DESCRIPTION
The one-sided crossing minimization problem (OCM) asks for an order of one layer of a two-layer bipartite graph, 
where the other layer of nodes is fixed, such that the number of edge crossings is minimized. 
Our heiCross solver employs a branch and reduce approach specifically designed for the heuristic track of the challenge. 

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

AUTHORS
Marlon Dittes @MarlonDittes
Shai Peretz @Shaip161
Antonie Wagner @Tonyloni

Co-Authors:
Alvaro Garmendia @alvarogarmen
Tomer Haham @TomerHaham
Ernestine Großman
Henning Woydt
Henrik Reinstädtler
Adil Chhabra