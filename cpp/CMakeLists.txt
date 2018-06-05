cmake_minimum_required(VERSION 3.10)

project(Hermite)
set (PROJECT_VERSION_MAJOR 0)
set (PROJECT_VERSION_MINOR 0)
set (CMAKE_VERBOSE_MAKEFILE OFF)
set (CMAKE_BUILD_TYPE release)
set (CMAKE_CXX_FLAGS_DEBUG "-DDEBUG")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -Wall -Wextra -Wpedantic") # Necessary to repeat
set (CMAKE_CXX_LINKER_FLAGS  "${CMAKE_CXX_LINKER_FLAGS} -Wl,--no-undefined")
set(CMAKE_CXX_STANDARD 14)

# Find Python
find_package(PythonLibs 3 REQUIRED)
include_directories(SYSTEM ${PYTHON_INCLUDE_DIRS})

# Find Boost
# set (BOOST_INCLUDEDIR ~/.guix-profile/include)
# set (BOOST_LIBRARYDIR ~/.guix-profile/lib)
find_package(Boost COMPONENTS python3 numpy3 REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

# Find sources and add include directory of sources
file(GLOB_RECURSE Sources "src/hermite/*.cpp")
include_directories("src")

# Build our library
add_library(hermite_cpp SHARED ${Sources})
target_link_libraries(hermite_cpp ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${CMAKE_DL_LIBS})

# Don't prepend library name with lib
set_target_properties(hermite_cpp PROPERTIES PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/..)

# Tests
include(CTest)
file(GLOB testfiles "src/tests/*.cpp")
foreach(testfile ${testfiles})
    get_filename_component(testname ${testfile} NAME_WE)
    add_executable(${testname}.exe ${testfile})
    target_link_libraries (${testname}.exe hermite_cpp)
    # set_target_properties(${testname}.exe PROPERTIES
    #     RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY})
    add_test(${testname} ${testname}.exe)
    set_tests_properties (${testname} PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
endforeach(testfile ${testfiles})

# Generate documentation
find_package(Doxygen)
set(DOXYGEN_CALL_GRAPH YES)
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_GENERATE_LATEX YES)
set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/doc)
set(DOXYGEN_PROJECT_NAME ${PROJECT_NAME}-${PROJECT_VERSION_MAJOR}-${PROJECT_VERSION_MINOR})
set(DOXYGEN_RECURSIVE YES)
set(DOXYGEN_USE_MATHJAX YES)
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
doxygen_add_docs("doc")
