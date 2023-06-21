cmake_minimum_required(VERSION 3.23)
project(MO10)

set(CMAKE_CXX_STANDARD 23)

add_executable(MO10 main.cpp nonlinearSolvers.h nonlinearSolvers.cpp nonlinearSolvers.h)
