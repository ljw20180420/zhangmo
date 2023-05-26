cmake_minimum_required(VERSION 3.22)
project(LPM)

add_executable(LPM main.cpp)
# set(CMAKE_BUILD_TYPE "Debug")
# target_compile_options(LPM PRIVATE -O0 -W -Wall -std=c++20)
set(CMAKE_BUILD_TYPE "Release")
target_compile_options(LPM PRIVATE -O3 -W -Wall -std=c++20)

install(TARGETS LPM RUNTIME DESTINATION bin)
