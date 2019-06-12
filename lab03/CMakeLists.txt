cmake_minimum_required(VERSION 2.8.9)
project (program)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
<<<<<<< HEAD
set(CMAKE_CXX_FLAGS "-Wall -g")
=======
set(CMAKE_CXX_FLAGS "-Wall -g -Ofast")
>>>>>>> 7d6ab7ed33b68910595c18ac5ec1d0589ce930fe

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME})
