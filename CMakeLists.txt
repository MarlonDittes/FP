cmake_minimum_required(VERSION 3.13)
project(test VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Release)

# Enable debugging information
#set(CMAKE_BUILD_TYPE Debug)

# Enable profiling
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

#add_library(my_lib STATIC graph.cpp io.cpp stack.cpp reductions.cpp performance.cpp LPSolver.cpp)
#add_library(my_lib STATIC graph.cpp io.cpp stack.cpp reductions.cpp performance.cpp)

add_executable(PACE2024Exact
    main.cpp
    graph.cpp
    io.cpp
    stack.cpp
    reductions.cpp
    performance.cpp
    src_henning/src/misc.cpp
    src_henning/src/edge.cpp
    src_henning/src/macros.h
    src_henning/src/definitions.h
    src_henning/src/misc.h
    src_henning/src/edge.h
    src_henning/src/graph_hen.h
    src_henning/src/solver_bf.h
    src_henning/src/perfect_pattern.h
    src_henning/src/exhaustive_solver.h
    src_henning/src/candidate_manager.h
    src_henning/src/cross_matrix.h
    src_henning/src/partitioner.h
    src_henning/src/translation_table.h
    src_henning/src/twin_reducer.h
    src_henning/src/solver.h
    src_henning/src/front_back_reducer.h
    src_henning/src/useless_reducer.h
    src_henning/src/one_vertex_partitioner.h
    src_henning/src/one_vertex_partition_manager.h
)

#add_executable(main main.cpp)
#target_link_libraries(main PRIVATE my_lib)

