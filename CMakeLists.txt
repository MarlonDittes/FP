cmake_minimum_required(VERSION 3.13)
project(test VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE RELEASE)


add_library(my_lib STATIC argtable3.c graph.cpp io.cpp)

add_executable(main main.cpp)
add_executable(mainarg mainarg.cpp)
target_link_libraries(main PRIVATE my_lib)
target_link_libraries(mainarg PRIVATE my_lib)

