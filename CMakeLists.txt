cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(bd_rate LANGUAGES CXX C)

include(FetchContent)

# setup Conan package finding.
set(CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})

# language settings.
set(CMAKE_CXX_STANDARD 17)  # C++17
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS "YES")

# alglib library
FetchContent_Declare(
    alglib_src
    URL https://www.alglib.net/translator/re/alglib-4.01.0.cpp.gpl.tgz
)
FetchContent_MakeAvailable(alglib_src)

message("alglib_src_SOURCE_DIR: ${alglib_src_SOURCE_DIR}")
message("alglib_src_BINARY_DIR: ${alglib_src_BINARY_DIR}")
file(GLOB alglib_SRC_FILES
    ${alglib_src_SOURCE_DIR}/src/*.cpp
    ${alglib_src_SOURCE_DIR}/src/*.h
)

add_library(alglib STATIC ${alglib_SRC_FILES})
target_compile_options(
    alglib PRIVATE
    -O3
    -DAE_CPU=AE_INTEL
    -mavx2
    -mfma
)

# main application.
find_package(cxxopts)

set(${PROJECT_NAME}_SRC_FILES
    bd_rate.cpp
)
add_executable(${PROJECT_NAME} ${bd_rate_SRC_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE alglib)
target_include_directories(${PROJECT_NAME} PUBLIC ${cxxopts_INCLUDE_DIR} ${alglib_src_SOURCE_DIR}/src/)
target_compile_options(
    ${PROJECT_NAME} PRIVATE
    -O3
    -Wall
    -Wextra
)

