cmake_minimum_required(VERSION 3.16)

project (LaticeBoltzmannSolver2D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SOURCES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INCLUDES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(SOURCES
    "${SOURCES_DIR}/BoundingBox.cpp"
    "${SOURCES_DIR}/main.cpp"
)

#TODO: Fix include filepaths
#TODO: Refer to old Makefile and add compiler flags
#TODO: Execute for all 3 parameter cases and measure runtimes

include_directories("${INCLUDES_DIR}")
include_directories("${SOURCES_DIR}")
add_executable(lbm ${SOURCES})


