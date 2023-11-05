cmake_minimum_required(VERSION 3.0.0)
project(practical VERSION 0.1.0)
set(CMAKE_BUILD_TYPE Release)

add_library(my_lib STATIC argtable3.c functions.cpp)

add_executable(main main.cpp)

target_link_libraries(main my_lib)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
