cmake_minimum_required(VERSION 3.22.1)

project (proj LANGUAGES CXX)

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


# # ##########
# Flint
include(FindPkgConfig)
pkg_check_modules(FLINT flint)

# include(FindPackageHandleStandardArgs)
link_libraries(${pkgcfg_lib_FLINT_flint} ${pkgcfg_lib_FLINT_gmp} ${pkgcfg_lib_FLINT_flint_mpfr})

include_directories("/usr/local/include/flint")


#
file(GLOB_RECURSE COMMON_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/common/*)
add_library(libcomp STATIC ${COMMON_SRC})
target_include_directories(libcomp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/common)
target_link_libraries(libcomp PRIVATE spu-static-lib)

add_executable(proj ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
target_link_libraries(proj PRIVATE libcomp)
target_link_libraries(proj PRIVATE spu-static-lib)
target_include_directories(proj PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/common)




# # Check for StarPU
# if (NOT STARPU_FIND_VERSION)
#   set (STARPU_FIND_VERSION "1.3")
#   message(STATUS "Checking for StarPU-${STARPU_FIND_VERSION} [default]")
# else()
#   message(STATUS "Checking for StarPU-${STARPU_FIND_VERSION}")
# endif()

