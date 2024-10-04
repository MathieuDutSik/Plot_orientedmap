cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Set the C++ standard to C++11 (or higher if needed)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 (make sure you have Eigen installed and CMake can find it)
find_package(Eigen3 3.3 REQUIRED)

# Find Boost (adjust components if you need specific ones like filesystem, system, etc.)
find_package(Boost REQUIRED COMPONENTS system filesystem)

# Add include directories (adjust paths to your actual include locations)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_basic)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_number)
include_directories(${CMAKE_SOURCE_DIR}/permutalib/src)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Add source files (replace with your actual source files if needed)
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Add the executable
add_executable(CombPlaneToSVG CombPlaneToSVG.cpp)

# Eigen is header-only, no need to link it explicitly, but its headers are included.
