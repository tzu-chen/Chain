cmake_minimum_required(VERSION 3.17)
project(ChainDMRG)

set(CMAKE_CXX_STANDARD 17)

set(ITENSOR_DIR /home/tzuchen/Software/itensor/)

include_directories(${ITENSOR_DIR}
        ${ITENSOR_DIR}itensor/
        ${ITENSOR_DIR}itensor/detail
        ${ITENSOR_DIR}itensor/itdata
        ${ITENSOR_DIR}itensor/mps
        ${ITENSOR_DIR}itensor/mps/lattice/
        ${ITENSOR_DIR}itensor/mps/sites/
        ${ITENSOR_DIR}itensor/tensor
        ${ITENSOR_DIR}itensor/util
        )

add_library(itensor STATIC IMPORTED)
set_property(TARGET itensor PROPERTY IMPORTED_LOCATION ${ITENSOR_DIR}lib/libitensor.a)
set(ITENSOR_FLAGS " -lpthread -I/usr/include/openblas -lopenblas -DPLATFORM_openblas -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE -L/home/tzuchen/Software/itensor/itensor -litensor")
set(ITENSOR_LINK_FLAGS "-DPLATFORM_openblas -fopenmp -L/home/tzuchen/Software/itensor/lib -litensor -lpthread -I/usr/lib64/openblas -lopenblas")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ITENSOR_FLAGS}")

add_executable(Chain dmrg.cpp f_data.cpp golden_site.cpp haagerup_q_site.cpp)
target_link_libraries(Chain itensor "${ITENSOR_LINK_FLAGS}")
