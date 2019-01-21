cmake_minimum_required(VERSION 3.5)

project(Hermite)
set (PROJECT_VERSION_MAJOR 0)
set (PROJECT_VERSION_MINOR 2)
set (CMAKE_VERBOSE_MAKEFILE OFF)
set (CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_RELEASE} -DDEBUG")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -Wall -Wextra -Wpedantic") # Necessary to repeat
set (CMAKE_CXX_LINKER_FLAGS  "${CMAKE_CXX_LINKER_FLAGS} -Wl,--no-undefined")
set (CMAKE_CXX_STANDARD 14)

# Default build type
#
# Use command-line argument to change:
#   cmake -D CMAKE_BUILD_TYPE=debug ...
if (CMAKE_BUILD_TYPE STREQUAL "")
    set (CMAKE_BUILD_TYPE release)
endif (CMAKE_BUILD_TYPE STREQUAL "")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Default library output directory
if (NOT DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    message(STATUS "Setting default library output directory...")
    set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/..)
endif (NOT DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
message(STATUS "Library output directory: '${CMAKE_LIBRARY_OUTPUT_DIRECTORY}'")

# Alternative to target only a specific library
# set_target_properties(hermite_cpp PROPERTIES PREFIX ""
#     LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/..)

# Find Python
find_package(PythonLibs 3 REQUIRED)
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})
message(STATUS "Python include directory: ${PYTHON_INCLUDE_DIRS}")

# Link Boost statically
set (BOOST_LIBRARYDIR ${CMAKE_HOME_DIRECTORY}/boost/lib)

if (NOT EXISTS /usr/include/boost/python/numpy)
    set (BOOST_INCLUDEDIR /usr/local/lib/R/site-library/BH/include)
endif (NOT EXISTS /usr/include/boost/python/numpy)

find_package(Boost COMPONENTS python3 numpy3)

# Fix for archlinux
if (NOT Boost_FOUND)
    find_package(Boost COMPONENTS python37 numpy37 REQUIRED)
endif (NOT Boost_FOUND)

message(STATUS "Boost include directory: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost library directory: ${Boost_LIBRARY_DIRS}")
# include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# Find sources and add include directory of sources
file(GLOB_RECURSE Sources "src/hermite/*.cpp")
include_directories("src")

# Build our library
add_library(hermite_cpp SHARED ${Sources})
target_link_libraries(hermite_cpp ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${CMAKE_DL_LIBS})

# Remove prefix
set_target_properties(hermite_cpp PROPERTIES PREFIX "")

# Tests
include(CTest)
file(GLOB testfiles "src/tests/*.cpp")
foreach(testfile ${testfiles})
    get_filename_component(testname ${testfile} NAME_WE)
    add_executable(${testname}.exe ${testfile})
    target_link_libraries (${testname}.exe hermite_cpp)
    add_test(${testname} ${testname}.exe)
    set_tests_properties (${testname} PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
endforeach(testfile ${testfiles})

# Files used to test ideas and library functions
file(GLOB scratchfiles "src/scratchpad/*.cpp")
foreach(scratchfile ${scratchfiles})
    get_filename_component(name ${scratchfile} NAME_WE)
    add_executable(${name}.exe ${scratchfile})
    target_link_libraries (${name}.exe hermite_cpp)
endforeach(scratchfile ${scratchfiles})
