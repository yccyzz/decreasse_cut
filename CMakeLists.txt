# CMake最低版本要求
cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER "gcc.exe")
set(CMAKE_CXX_COMPILER "g++.exe")

# 设置C++标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 项目名称和版本
project(FPGAReader)


# 设置编译器标志
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")


# 源文件
set(SOURCES
    fpga_reader.cpp
        main.cpp
        fpga_reader.h
)

# 创建可执行文件
add_executable(fpga_reader ${SOURCES})

