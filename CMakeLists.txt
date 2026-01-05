cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(CephSim)

# Search for SimGrid
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(SimGrid 3.27 REQUIRED) # This template requires SimGrid v3.27 because we use templated get()

include_directories(${SimGrid_INCLUDE_DIR})

# add_executable(ping-pong ping-pong.cpp)
# target_link_libraries(ping-pong ${SimGrid_LIBRARY})

# add_executable(master-workers master-workers.cpp)
# target_link_libraries(master-workers ${SimGrid_LIBRARY})

# my stuff
add_executable(ceph-sim
    src/ceph-sim.cpp
    src/SimContext.cpp
    src/CephCommon.cpp
    src/OsdActor.cpp
    src/MonActor.cpp
    src/CephActor.cpp
    src/ClientActor.cpp
    src/MetricMonitor.cpp
)
target_include_directories(ceph-sim PUBLIC src)
target_link_libraries(ceph-sim PRIVATE ${SimGrid_LIBRARY})

# Integrate dmclock
include(FetchDmclock)
target_link_libraries(ceph-sim PRIVATE dmclock::dmclock)
