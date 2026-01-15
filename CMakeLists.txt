cmake_minimum_required(VERSION 3.10...3.25)

project(negative_triangle_counting)

# Set C++ version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

#Sanitizer options
# option(ENABLE_SANITIZERS "Build with Address/Undefined sanitizers" ON)

if (ENABLE_SANITIZERS)
    message(STATUS "Sanitizers enabled")
    set(SANITIZER_FLAGS "-fsanitize=address,undefined -fno-omit-frame-pointer -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZER_FLAGS}")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${SANITIZER_FLAGS}")
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif ()

# Enable OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

find_package(Boost REQUIRED)

include_directories(include)
include_directories(${CMAKE_SOURCE_DIR}/external)

add_executable(benchmark
        src/benchmark.cpp
        src/counting.cpp
        src/graph_generation.cpp
        src/utils.cpp
        src/smooth_sensitivity.cpp
        src/global_sensitivity.cpp
        src/graph_statistics.cpp
        src/smooth_sensitivity.cpp
        src/smooth_sensitivity.cpp
        src/smooth_sensitivity_unbiased.cpp
        include/smooth_sensitivity_unbiased.h
        src/distance_retrieval.cpp
        include/distance_retrieval.h
        src/smooth_sens_utils.cpp
        include/smooth_sens_utils.h
)

target_link_libraries(benchmark
        Boost::headers)

# Enable testing support
enable_testing()

#add_executable(treap_test test/treap_test.cpp src/counting.cpp src/graph_generation.cpp src/utils.cpp src/smooth_sensitivity.cpp src/global_sensitivity.cpp
#        src/smooth_sensitivity.cpp
#        src/smooth_sensitivity.cpp)
#target_link_libraries(treap_test
#        Boost::headers
#)
#
#add_executable(counting_test test/counting_test.cpp src/counting.cpp src/graph_generation.cpp src/utils.cpp src/smooth_sensitivity.cpp src/global_sensitivity.cpp src/graph_statistics.cpp
#        src/smooth_sensitivity.cpp
#        src/smooth_sensitivity.cpp)
#target_link_libraries(counting_test
#        Boost::headers
#)
#
#add_executable(assignment_test test/assignment_test.cpp src/counting.cpp src/graph_generation.cpp src/utils.cpp src/smooth_sensitivity.cpp src/global_sensitivity.cpp src/graph_statistics.cpp
#        src/smooth_sensitivity.cpp
#        src/smooth_sensitivity.cpp)
#target_link_libraries(assignment_test
#        Boost::headers
#)
#
#add_executable(sensitivity_test test/sensitivity_test.cpp src/counting.cpp src/graph_generation.cpp src/utils.cpp src/smooth_sensitivity.cpp src/global_sensitivity.cpp src/graph_statistics.cpp
#        src/smooth_sensitivity.cpp
#        src/smooth_sensitivity.cpp)
#target_link_libraries(sensitivity_test
#        Boost::headers
#)


add_executable(distance_test
        test/distance_test.cpp
        src/distance_retrieval.cpp
        test/smooth_sens_unb_test.cpp
        test/smooth_sens_unb_test.h
        test/graph_presets.h
)

target_link_libraries(distance_test
        Boost::headers
)

if (OpenMP_CXX_FOUND)
    #    target_link_libraries(treap_test OpenMP::OpenMP_CXX)
    #    target_link_libraries(assignment_test OpenMP::OpenMP_CXX)
    #    target_link_libraries(sensitivity_test OpenMP::OpenMP_CXX)
    #    target_link_libraries(counting_test OpenMP::OpenMP_CXX)
    #    target_link_libraries(distribution_test OpenMP::OpenMP_CXX)
    target_link_libraries(benchmark OpenMP::OpenMP_CXX)
    target_link_libraries(distance_test OpenMP::OpenMP_CXX)
endif ()

#add_test(NAME TreapTest COMMAND treap_test)
#add_test(NAME SensitivityTest COMMAND sensitivity_test)
#add_test(NAME AssignmentTest COMMAND assignment_test)
#add_test(NAME CountingTest COMMAND counting_test)
#add_test(NAME DistributionTest COMMAND distribution_test)
