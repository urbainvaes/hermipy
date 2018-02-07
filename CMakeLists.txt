cmake_minimum_required(VERSION 3.10)

project(mckean-vlasov)
set (mckean-vlasov_VERSION_MAJOR 1)
set (mckean-vlasov_VERSION_MINOR 0)
 
# Find necessary packages
find_package(PythonLibs 3 REQUIRED)
find_package(Boost COMPONENTS python3 REQUIRED)
include_directories(${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})
 
# Build our library
add_library(hermite SHARED integrator.cpp)
target_link_libraries(hermite ${Boost_LIBRARIES})

# Don't prepend library name with lib
set_target_properties(hermite PROPERTIES
    PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}
    )
