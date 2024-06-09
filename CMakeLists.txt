cmake_minimum_required(VERSION 3.13)
project(test VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)

# Need this for optil static binary submissions
set(CMAKE_EXE_LINKER_FLAGS "-static")

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

#add_library(my_lib STATIC src/graph.cpp src/io.cpp src/stack.cpp src/reductions.cpp src/branchandreduce.cpp)
#add_executable(heiCross main.cpp)
#target_link_libraries(heiCross PRIVATE my_lib)

add_executable(heiCross
    main.cpp
    src/graph.cpp
    src/io.cpp
    src/stack.cpp
    src/reductions.cpp
    src/branchandreduce.cpp
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