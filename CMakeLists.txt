cmake_minimum_required(VERSION 3.17)
project(ray_tracing)
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3")
include_directories(.)

add_executable(ray_tracing
        ray_tracing.cpp
        geometry.h
        Objects.hpp
        Objects.cpp Time.hpp Time.cpp  Material.cpp Material.hpp Lights.hpp Lights.cpp)
