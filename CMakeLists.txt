cmake_minimum_required(VERSION 4.0)
project(Spark)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

set(spark_SOURCE
    src/cmdline.cc
    src/sparse_tree.cc
    src/Spark.cc
    src/trace_arrow.cc 
)
set(vienna_SOURCE
    src/ViennaRNA/params/default.cc
    src/ViennaRNA/params/params.cc
    src/ViennaRNA/params/io.cc
    src/ViennaRNA/model.cc 
    src/ViennaRNA/utils.cc
)

add_library(RNA STATIC ${vienna_SOURCE})

target_include_directories(RNA PRIVATE .)

add_executable(Spark ${spark_SOURCE})

target_link_libraries(Spark PRIVATE RNA)

include_directories(src)
