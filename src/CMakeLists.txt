cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(OctreeBspline)

find_package(Eigen3 3.3.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(VTK 5.8 REQUIRED)
include(${VTK_USE_FILE})

file(GLOB OctreeBspline_SRC
    "*.h"
    "*.c"
    "*.cpp"
    "*.inl"
)

add_executable(OctreeBspline ${OctreeBspline_SRC})

target_link_libraries(OctreeBspline ${VTK_LIBRARIES})