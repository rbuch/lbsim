cmake_minimum_required(VERSION 3.10)

project(lbsim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

find_library(METIS metis)

add_executable(lbsim src/main.cpp)
target_link_libraries(lbsim PRIVATE ${METIS})

add_executable(constrain src/constrain.cpp)
target_link_libraries(constrain PRIVATE ${METIS})
