cmake_minimum_required(VERSION 3.8...3.28)
project(KarstNSim VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Engine sources. main.cpp (the interactive/CLI front end) is not part of the library.
set(KARSTNSIM_CORE_SOURCES
    src/karstic_network.cpp
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
    src/job_context.cpp
    src/library.cpp
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
    include/KarstNSim/errors.h
    include/KarstNSim/job_context.h
    include/KarstNSim/library.h
)

# Reusable native backend. Static by default (BUILD_SHARED_LIBS selects shared);
# position-independent so it can be linked into shared modules.
add_library(karstnsim_core ${KARSTNSIM_CORE_SOURCES} ${KARSTNSIM_INCLUDES})
target_include_directories(karstnsim_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_features(karstnsim_core PUBLIC cxx_std_17)
set_target_properties(karstnsim_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Keep intra-library inlining identical to the non-PIC executable build.
    target_compile_options(karstnsim_core PRIVATE -fno-semantic-interposition)
endif()

# Legacy command-line executable, unchanged behavior and parameters.
add_executable(karstnsim src/main.cpp)
target_link_libraries(karstnsim PRIVATE karstnsim_core)

include(CTest)
option(KARSTNSIM_SANITIZERS "Enable address and undefined-behavior sanitizers" OFF)
set(KARSTNSIM_TARGETS karstnsim_core karstnsim)
if(BUILD_TESTING)
    add_executable(packed_adjacency_test tests/packed_adjacency_test.cpp)
    target_include_directories(packed_adjacency_test PRIVATE include)
    add_test(NAME packed_adjacency COMMAND packed_adjacency_test)
    add_executable(cost_graph_test tests/cost_graph_test.cpp)
    target_link_libraries(cost_graph_test PRIVATE karstnsim_core)
    add_test(NAME cost_graph COMMAND cost_graph_test)
    add_executable(compact_sampler_grid_test tests/compact_sampler_grid_test.cpp)
    target_include_directories(compact_sampler_grid_test PRIVATE include)
    add_test(NAME compact_sampler_grid COMMAND compact_sampler_grid_test)
    find_package(Threads REQUIRED)
    add_executable(native_library_test tests/native_library_test.cpp)
    target_link_libraries(native_library_test PRIVATE karstnsim_core Threads::Threads)
    add_test(NAME native_library COMMAND native_library_test)
    set_tests_properties(native_library PROPERTIES TIMEOUT 600)
    list(APPEND KARSTNSIM_TARGETS packed_adjacency_test cost_graph_test compact_sampler_grid_test native_library_test)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
        add_test(NAME regression_exports COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/run_regression.py
            --binary $<TARGET_FILE:karstnsim>
            --out ${CMAKE_CURRENT_BINARY_DIR}/regression-results
            --golden ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/golden/expected_optimized.json)
        set_tests_properties(regression_exports PROPERTIES TIMEOUT 600
            FIXTURES_SETUP legacy_regression_exports)
        # Reproduces the legacy exports just checked above through run_simulation_memory.
        add_executable(native_parity_test tests/native_parity_test.cpp)
        target_link_libraries(native_parity_test PRIVATE karstnsim_core)
        list(APPEND KARSTNSIM_TARGETS native_parity_test)
        add_test(NAME native_parity COMMAND native_parity_test
            ${CMAKE_CURRENT_BINARY_DIR}/regression-results)
        set_tests_properties(native_parity PROPERTIES TIMEOUT 600
            FIXTURES_REQUIRED legacy_regression_exports)
    endif()
endif()

foreach(target IN LISTS KARSTNSIM_TARGETS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
        if(KARSTNSIM_SANITIZERS)
            target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer)
            target_link_libraries(${target} PRIVATE $<BUILD_INTERFACE:-fsanitize=address,undefined>)
        endif()
    endif()
endforeach()
