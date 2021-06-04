## If using older version, change to your CMake version and see if it works.
cmake_minimum_required(VERSION 2.8)
project(ChainDMRG)

## Tells the compiler to use c++17 standard.
if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17 -fopenmp")
    else()
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")
    endif ()
else ()
    set (CMAKE_CXX_STANDARD 17)
endif ()

## Default project structure. Shouldn't have to modify
include_directories( ./include ./src extern/itensor/include ./extern)
add_library(itensor STATIC IMPORTED)

## Your blas header and library location. Depends on OS and installation
set(BLAS_LAPACK_FLAG "-framework Accelerate")
set(BLAS_LAPACK_INCLUDEFLAGS "")

set(ITENSOR_FLAGS "${BLAS_LAPACK_INCLUDEFLAGS} -I${CMAKE_SOURCE_DIR}/extern/itensor/include/itensor")
set(ITENSOR_LINK_FLAGS "-L${CMAKE_SOURCE_DIR}/extern/itensor/lib -litensor ${BLAS_LAPACK_FLAG}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ITENSOR_FLAGS} -DPLATFORM_macos")

add_executable(Chain src/dmrg.cpp src/f_data.cpp src/golden_site.cpp src/haagerup_q_site.cpp src/chain.cpp src/haagerup_site.cpp src/utility.cpp)

set_property(TARGET itensor PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/extern/itensor/lib/libitensor.a)
target_link_libraries(Chain itensor  "${ITENSOR_LINK_FLAGS}")
