cmake_minimum_required(VERSION 3.8...3.28)
project(KarstNSim VERSION 1.0 LANGUAGES CXX)



set(CMAKE_CXX_STANDARD 14 CACHE STRING "The C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(KARSTNSIM_SOURCES src/KarsticNetwork.cpp src/main.cpp src/surface_sampling.cpp src/cost_graph.cpp src/randomgenerator.cpp  src/write_files.cpp src/graph_operations.cpp  src/read_files.cpp
src/karstic_skeleton.cpp  src/run_code.cpp src/parse_inputs.cpp
)
set(KARSTNSIM_INCLUDES include/KarstNSim/KarsticNetwork.h include/KarstNSim/surface_sampling.h include/KarstNSim/basics.h include/KarstNSim/randomgenerator.h  include/KarstNSim/vec.h
include/KarstNSim/geology.h         include/KarstNSim/read_files.h       include/KarstNSim/write_files.h
include/KarstNSim/graph.h           include/KarstNSim/run_code.h include/KarstNSim/parse_inputs.h
)

include_directories("include")

add_executable(karstnsim ${KARSTNSIM_SOURCES} ${KARSTNSIM_INCLUDES})
