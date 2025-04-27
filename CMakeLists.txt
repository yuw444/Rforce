cmake_minimum_required(VERSION 3.16.0)

set(CMAKE_C_COMPILER /home/yu89975/local/bin/gcc)
set(CMAKE_CXX_COMPILER /home/yu89975/local/bin/g++)

project(Rforce VERSION 1.0)

set(CMAKE_C_STANDARD 99)
set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Release)

find_package(OpenMP REQUIRED)

include_directories(./criterion/include /home/yu89975/local/include /hpc/apps/R/4.3.3/lib64/R/include)
link_directories(./criterion/lib /home/yu89975/local/lib /hpc/apps/R/4.3.3/lib64/R/lib)

set(SOURCES
    ./src/argparse.c
    ./src/forest.c 
    ./src/tree.c 
    ./src/split.c
    ./src/utils.c 
    ./src/surv.c
    ./src/gee.c
    ./src/pvalue.c
    ./src/main.c
)

add_executable(Rforce ${SOURCES})
target_link_libraries(Rforce PRIVATE OpenMP::OpenMP_C gsl gslcblas m g)


