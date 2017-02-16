cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(OctreeBspline)

file(GLOB OctreeBspline_SRC
    "*.h"
    "*.cpp"
    "*.inl"
)

add_executable(OctreeBspline ${OctreeBspline_SRC})