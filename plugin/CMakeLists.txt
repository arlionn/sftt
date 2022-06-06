cmake_minimum_required(VERSION 3.14)
project(stplugin)

set(CMAKE_CXX_STANDARD 11)
set(BOOST_ROOT "/usr/local/include/boost")
include_directories(/usr/local/Cellar/boost/1.71.0/include)

add_executable(main main.cpp )
add_executable(testboost main.cpp owen.hpp owen.cpp)