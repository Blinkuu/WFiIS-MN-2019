cmake_minimum_required(VERSION 2.8.9)
project (program)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -g")

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME})
