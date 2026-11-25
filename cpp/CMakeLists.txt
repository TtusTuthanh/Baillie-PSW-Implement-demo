cmake_minimum_required(VERSION 3.10)
project(bpsw_cpp VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# __int128 is a GCC/Clang extension; MSVC doesn't support it
if(MSVC)
    message(FATAL_ERROR "MSVC does not support __int128. Use GCC or Clang.")
endif()

add_library(bpsw
    bpsw.cpp
)

target_compile_options(bpsw PRIVATE -Wall -Wextra -O2)

add_executable(bpsw_test test_bpsw.cpp)
target_link_libraries(bpsw_test PRIVATE bpsw)

include(CTest)
enable_testing()
add_test(NAME bpsw_tests COMMAND bpsw_test)
