cmake_minimum_required (VERSION 3.14)

include(FetchContent)
include(CheckLanguage)

project (rref LANGUAGES CXX C)

include(GNUInstallDirs)

set( CMAKE_CXX_FLAGS "-std=c++17" )
add_compile_options( -g )
add_compile_options( -Wfatal-errors )

add_library(rref SHARED rref.cpp )
set_target_properties( rref PROPERTIES LIBRARY_OUTPUT_DIRECTORY "." )


add_executable( test test.cpp )
target_link_libraries( test rref )
