cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(SimGridTemplateProject) # TODO: give a real name to your project here

# Search for SimGrid
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(SimGrid 3.27 REQUIRED) # This template requires SimGrid v3.27 because we use templated get()

include_directories(${SimGrid_INCLUDE_DIR})

add_executable(ping-pong ping-pong.cpp)
target_link_libraries(ping-pong ${SimGrid_LIBRARY})

add_executable(master-workers master-workers.cpp)
target_link_libraries(master-workers ${SimGrid_LIBRARY})

