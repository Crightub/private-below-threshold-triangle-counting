cmake_minimum_required(VERSION 3.10...3.25)

project(negative_triangle_counting)

# Set C++ version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(/usr/local/opt/boost/include)
include_directories(include)

find_package(nlohmann_json REQUIRED)


add_subdirectory(osqp-cpp)

add_executable(negative_triangle_counting negative_triangle_count.cpp src/basic_counting.cpp src/distribution.cpp src/graph_generation.cpp src/qp_counting.cpp src/utils.cpp)

# Link libraries
target_link_libraries(negative_triangle_counting
    osqp-cpp
)

target_link_libraries(negative_triangle_counting
    nlohmann_json::nlohmann_json
)

# OpenMP if needed
if(OpenMP_CXX_FOUND)
    target_link_libraries(negative_triangle_counting OpenMP::OpenMP_CXX)
endif()


add_executable(benchmark
    benchmark.cpp
    src/basic_counting.cpp
    src/distribution.cpp
    src/graph_generation.cpp
    src/qp_counting.cpp
    src/utils.cpp
)

target_link_libraries(benchmark
    osqp-cpp
)

target_link_libraries(benchmark nlohmann_json::nlohmann_json)

if(OpenMP_CXX_FOUND)
    target_link_libraries(benchmark OpenMP::OpenMP_CXX)
endif()
