cmake_minimum_required(VERSION 3.20)
project(acm2)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(AlgLib)

add_executable(acm2 main.cpp)

target_link_libraries(acm2 PRIVATE BigInt)
target_link_libraries(acm2 PRIVATE AlgLib)