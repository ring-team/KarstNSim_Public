cmake_minimum_required(VERSION 3.8...3.28)
project(KarstNSim VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(KARSTNSIM_SOURCES 
    src/karstic_network.cpp 	
    src/main.cpp 
    src/surface_sampling.cpp 
    src/cost_graph.cpp 
    src/geostats.cpp 
    src/ghost_rocks.cpp 
    src/randomgenerator.cpp  
    src/write_files.cpp 
    src/simplex_noise.cpp 
    src/graph_operations.cpp  
    src/read_files.cpp
    src/karstic_skeleton.cpp  
    src/run_code.cpp
    src/parse_inputs.cpp
)

set(KARSTNSIM_INCLUDES 
    include/KarstNSim/karstic_network.h 
    include/KarstNSim/surface_sampling.h 
    include/KarstNSim/basics.h 
    include/KarstNSim/randomgenerator.h  
    include/KarstNSim/vec.h
    include/KarstNSim/geology.h         
    include/KarstNSim/geostats.h        
    include/KarstNSim/ghost_rocks.h    
    include/KarstNSim/simplex_noise.h    
    include/KarstNSim/read_files.h       
    include/KarstNSim/write_files.h
    include/KarstNSim/graph.h           
    include/KarstNSim/run_code.h 
    include/KarstNSim/parse_inputs.h 
    include/KarstNSim/nanoflann.hpp
    include/KarstNSim/models/results.h
)

include_directories("include")

add_executable(karstnsim ${KARSTNSIM_SOURCES} ${KARSTNSIM_INCLUDES})

include(CTest)
option(KARSTNSIM_SANITIZERS "Enable address and undefined-behavior sanitizers" OFF)
set(KARSTNSIM_TARGETS karstnsim)
if(BUILD_TESTING)
    add_executable(packed_adjacency_test tests/packed_adjacency_test.cpp)
    add_test(NAME packed_adjacency COMMAND packed_adjacency_test)
    add_executable(cost_graph_test tests/cost_graph_test.cpp src/cost_graph.cpp)
    add_test(NAME cost_graph COMMAND cost_graph_test)
    add_executable(compact_sampler_grid_test tests/compact_sampler_grid_test.cpp)
    add_test(NAME compact_sampler_grid COMMAND compact_sampler_grid_test)
    list(APPEND KARSTNSIM_TARGETS packed_adjacency_test cost_graph_test compact_sampler_grid_test)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
        add_test(NAME regression_exports COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/run_regression.py
            --binary $<TARGET_FILE:karstnsim>
            --out ${CMAKE_CURRENT_BINARY_DIR}/regression-results
            --golden ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/golden/expected_optimized.json)
        set_tests_properties(regression_exports PROPERTIES TIMEOUT 600)
    endif()
endif()

foreach(target IN LISTS KARSTNSIM_TARGETS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
        if(KARSTNSIM_SANITIZERS)
            target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
            target_link_libraries(${target} PRIVATE -fsanitize=address,undefined)
        endif()
    endif()
endforeach()
