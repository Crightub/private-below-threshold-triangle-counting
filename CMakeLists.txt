cmake_minimum_required(VERSION 3.10...3.25)

project(negative_triangle_counting)

# Set C++ version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Sanitizer options
option(ENABLE_SANITIZERS "Build with Address/Undefined sanitizers" ON)

if (ENABLE_SANITIZERS)
    message(STATUS "Sanitizers enabled")
    set(SANITIZER_FLAGS "-fsanitize=address,undefined -fno-omit-frame-pointer -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZER_FLAGS}")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${SANITIZER_FLAGS}")
endif ()


# Enable OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

include_directories(/usr/local/opt/boost/include)
include_directories(include)

find_package(nlohmann_json REQUIRED)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(Boost REQUIRED COMPONENTS log log_setup thread system)

# uncomment for build with quadratic optimizer
#add_subdirectory(osqp-cpp)
#add_executable(negative_triangle_counting negative_triangle_count.cpp src/basic_counting.cpp src/graph_generation.cpp src/qp_counting.cpp src/utils.cpp)

add_executable(negative_triangle_counting negative_triangle_count.cpp src/basic_counting.cpp src/graph_generation.cpp src/utils.cpp src/smooth_sensitivity.cpp src/global_sensitivity.cpp)

target_link_libraries(negative_triangle_counting
        nlohmann_json::nlohmann_json
        Eigen3::Eigen
        Boost::log
        Boost::log_setup
        Boost::thread
        Boost::system
)

add_executable(benchmark
        src/benchmark.cpp
        src/basic_counting.cpp
        src/graph_generation.cpp
        src/utils.cpp
        src/smooth_sensitivity.cpp
        src/global_sensitivity.cpp
)

target_link_libraries(benchmark
        nlohmann_json::nlohmann_json
        Eigen3::Eigen
        Boost::log
        Boost::log_setup
        Boost::thread
        Boost::system
)


# Enable testing support
enable_testing()

add_executable(triangle_test test/triangle_test.cpp src/basic_counting.cpp src/graph_generation.cpp src/utils.cpp src/smooth_sensitivity.cpp src/global_sensitivity.cpp)
target_link_libraries(triangle_test
        nlohmann_json::nlohmann_json
        Eigen3::Eigen
        Boost::log
        Boost::log_setup
        Boost::thread
        Boost::system
)

find_package(Python3 COMPONENTS Development NumPy REQUIRED)

add_executable(counting_test test/counting_test.cpp src/basic_counting.cpp src/graph_generation.cpp src/utils.cpp src/smooth_sensitivity.cpp src/global_sensitivity.cpp)

target_include_directories(counting_test PRIVATE
        ${Python3_INCLUDE_DIRS}
        ${Python3_NumPy_INCLUDE_DIRS}
)

target_link_libraries(counting_test
        nlohmann_json::nlohmann_json
        Eigen3::Eigen
        Boost::log
        Boost::log_setup
        Boost::thread
        Boost::system
        ${Python3_LIBRARIES}
)

# uncomment for build with quadratic optimizer
#target_link_libraries(negative_triangle_counting
#    osqp-cpp
#)

# OpenMP if needed
if (OpenMP_CXX_FOUND)
    target_link_libraries(negative_triangle_counting OpenMP::OpenMP_CXX)
    target_link_libraries(triangle_test OpenMP::OpenMP_CXX)
    target_link_libraries(counting_test OpenMP::OpenMP_CXX)
    target_link_libraries(benchmark OpenMP::OpenMP_CXX)
endif ()

add_test(NAME TriangleTest COMMAND triangle_test)
add_test(NAME CountingTest COMMAND counting_test)