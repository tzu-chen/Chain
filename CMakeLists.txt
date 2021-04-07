cmake_minimum_required(VERSION 3.17)
project(ChainDMRG)

set(CMAKE_CXX_STANDARD 17)

include_directories( ./include ./src ./libs/itensor/include )

add_library(itensor STATIC IMPORTED)

set(ITENSOR_FLAGS " -lpthread -I/usr/include/openblas -lopenblas -DPLATFORM_openblas -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE -L${CMAKE_SOURCE_DIR}/libs/itensor/include/itensor -litensor")
set(ITENSOR_LINK_FLAGS "-DPLATFORM_openblas -fopenmp -L${CMAKE_SOURCE_DIR}/libs/itensor/lib -litensor -lpthread -I/usr/lib64/openblas -lopenblas")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ITENSOR_FLAGS}")

add_executable(Chain src/dmrg.cpp src/f_data.cpp src/golden_site.cpp src/haagerup_q_site.cpp)

set_property(TARGET itensor PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/libs/itensor/lib/libitensor.a)
target_link_libraries(Chain itensor "${ITENSOR_LINK_FLAGS}")
