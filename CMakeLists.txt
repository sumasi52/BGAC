cmake_minimum_required(VERSION 3.13)
project(bdd_circuit_approximate)

set(CMAKE_CXX_STANDARD 14)

add_executable(bdd_circuit_approximate main.cpp miniBDD.h miniBDD.cpp miniBDD.inc ga.cpp ga.h)