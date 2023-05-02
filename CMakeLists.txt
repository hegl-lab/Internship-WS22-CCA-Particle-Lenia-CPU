cmake_minimum_required(VERSION 3.22)
project(particle_lenia_cpu)

set(CMAKE_CXX_STANDARD 14)

# Set this to the minimal version you want to support
find_package(raylib 3.0 QUIET) # Let CMake search for a raylib-config.cmake

# You could change the QUIET above to REQUIRED and remove this if() clause
# This part downloads raylib and builds it if it's not installed on your system
if (NOT raylib_FOUND) # If there's none, fetch and build raylib
    include(FetchContent)

    FetchContent_Declare(
            raylib
            URL https://github.com/raysan5/raylib/archive/master.tar.gz
    )

    FetchContent_GetProperties(raylib)
    if (NOT raylib_POPULATED) # Have we downloaded raylib yet?
        set(FETCHCONTENT_QUIET NO)
        FetchContent_Populate(raylib)

        set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # don't build the supplied examples

        # build raylib
        add_subdirectory(${raylib_SOURCE_DIR} ${raylib_BINARY_DIR})

    endif()

endif()


# TODO: introduce similar flags for different compilers
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("-- Compiler is GNU, adding compiler flags to improve performance")
    add_compile_options(-O3 -ffast-math -funsafe-math-optimizations -msse4.2)
endif ()

add_library(particle_lenia src/ParticleLeniaBucket.cpp src/ParticleLenia.cpp src/ParticleLeniaBucketThreaded.cpp)

add_executable(particle_lenia_cpu src/main.cpp)
add_executable(test_bucket src/test_bucket.cpp)
add_executable(gui_bucket src/gui_bucket.cpp)
add_executable(gui_simple src/gui_simple.cpp)

target_link_libraries(particle_lenia_cpu particle_lenia)
target_link_libraries(test_bucket particle_lenia)
target_link_libraries(gui_bucket particle_lenia raylib)
target_link_libraries(gui_simple particle_lenia raylib)
