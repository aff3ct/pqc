cmake_minimum_required(VERSION 3.22.1)

project (pqc LANGUAGES CXX)

# set CMAKE_INSTALL_BINDIR, CMAKE_INSTALL_LIBDIR, CMAKE_INSTALL_INCLUDEDIR and CMAKE_INSTALL_DATAROOTDIR variables
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## StreamPU runtime system

### Link with the "Threads library (required to link with StreamPU after)
set(CMAKE_THREAD_PREFER_PTHREAD ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/streampu/include/streampu.hpp")
    message(STATUS "Using sub-project StreamPU in ${CMAKE_CURRENT_SOURCE_DIR}/streampu/")
    option(SPU_COMPILE_STATIC_LIB "" ON)
    option(SPU_COMPILE_SHARED_LIB "" OFF)
    option(SPU_TESTS "" OFF)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/streampu/)
    set(PQC_SPU_LIB spu-static-lib)
else()
    message(STATUS "Looking for pre-installed StreamPU")
    find_package(cpptrace REQUIRED)
    find_package(StreamPU REQUIRED)
    set(PQC_SPU_LIB spu::spu-static-lib)
endif()

# ##########

# Flint Library
find_package(PkgConfig)
pkg_check_modules(FLINT REQUIRED flint)

# pqc lib
file(GLOB_RECURSE COMMON_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/*)
add_library(pqc-static-lib STATIC ${COMMON_SRC})

target_include_directories(pqc-static-lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/pqc}>)
target_link_libraries(pqc-static-lib ${PQC_SPU_LIB})
target_link_libraries(pqc-static-lib ${FLINT_STATIC_LIBRARIES})
target_link_options(pqc-static-lib PUBLIC ${FLINT_STATIC_LDFLAGS})

install(TARGETS pqc-static-lib
        EXPORT pqc-config
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/
        COMPONENT library)

install(EXPORT
        pqc-config
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/pqc"
        NAMESPACE pqc::
        COMPONENT library)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pqc
        COMPONENT headers
        FILES_MATCHING PATTERN "*.hpp")

# Examples
add_executable(pqc-example ${CMAKE_CURRENT_SOURCE_DIR}/examples/main.cpp)
target_link_libraries(pqc-example pqc-static-lib)

add_executable(pqc-bike ${CMAKE_CURRENT_SOURCE_DIR}/examples/main_bike.cpp)
target_link_libraries(pqc-bike PRIVATE pqc-static-lib)

add_executable(pqc-cm ${CMAKE_CURRENT_SOURCE_DIR}/examples/main_cm.cpp)
target_link_libraries(pqc-cm PRIVATE pqc-static-lib)

add_executable(pqc-hqc ${CMAKE_CURRENT_SOURCE_DIR}/examples/main_hqc.cpp)
target_link_libraries(pqc-hqc PRIVATE pqc-static-lib)

install(TARGETS pqc-bike pqc-cm pqc-hqc
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/pqc/examples"
    COMPONENT examples)

