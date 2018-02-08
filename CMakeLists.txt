cmake_minimum_required(VERSION 3.10)

project(Hermite)
set (PROJECT_VERSION_MAJOR 0)
set (PROJECT_VERSION_MINOR 0)

# Find necessary packages
include_directories("libhermite/src/")

# Find sources and add include directory of sources
file(GLOB_RECURSE Sources "*.cpp")
include_directories("src")

# Build our library
add_library(hermite SHARED ${Sources})
target_link_libraries(hermite ${Boost_LIBRARIES})

# Don't prepend library name with lib
set_target_properties(hermite PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}
    )

# Generate documentation
find_package(Doxygen)
set(DOXYGEN_CALL_GRAPH YES)
set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_GENERATE_LATEX YES)
set(DOXYGEN_GENERATE_LATEX YES)
set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/doc)
set(DOXYGEN_PROJECT_NAME ${PROJECT_NAME}-${PROJECT_VERSION_MAJOR}-${PROJECT_VERSION_MINOR})
set(DOXYGEN_RECURSIVE YES)
set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
doxygen_add_docs("doc")
