cmake_minimum_required(VERSION 3.22.1)

project (pqc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Link with the "Threads library (required to link with streampu after)
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# ##########
# streampu
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/streampu/include/streampu.hpp")
    option(SPU_COMPILE_STATIC_LIB "" ON)
    option(SPU_COMPILE_SHARED_LIB "" OFF)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/streampu/)

    message(STATUS "Sub-project found: streampu")
else()
    message(FATAL_ERROR "streampu can't be found, try to init the submodule with the following cmd:\n"
                        "$ git submodule update --init --recursive")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/streampu/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/streampu/lib/rang/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/streampu/lib/cpptrace/include)
link_libraries(spu-static-lib)

# ##########

# Flint Library
include(FindPkgConfig)
pkg_check_modules(FLINT flint)

link_libraries(${FLINT_LIBRARIES})
link_directories(${FLINT_LIBRARY_DIRS})
include_directories(${FLINT_INCLUDE_DIRS})

# pqc lib
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/common)

file(GLOB_RECURSE COMMON_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/common/*)
add_library(pqc-static-lib STATIC ${COMMON_SRC})

# Examples
add_executable(pqc-example ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
target_link_libraries(pqc-example PRIVATE pqc-static-lib)

add_executable(pqc-bike ${CMAKE_CURRENT_SOURCE_DIR}/src/main_bike.cpp)
target_link_libraries(pqc-bike PRIVATE pqc-static-lib)

add_executable(pqc-cm ${CMAKE_CURRENT_SOURCE_DIR}/src/main_cm.cpp)
target_link_libraries(pqc-cm PRIVATE pqc-static-lib)

add_executable(pqc-hqc ${CMAKE_CURRENT_SOURCE_DIR}/src/main_hqc.cpp)
target_link_libraries(pqc-hqc PRIVATE pqc-static-lib)

