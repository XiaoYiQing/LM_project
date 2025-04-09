cmake_minimum_required(VERSION 3.5)
project(MyFuncs)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


option(BUILD_SHARED_LIBS "Build shared libraries" OFF)


# Create our main executable.
add_executable( myapp src/main.cpp )



# ======================================================================= >>>>>
#   Present Project Additional Resource Inclusion
# ======================================================================= >>>>>

# Define the source directory and a destination directory as variables.
set(RES_XYQ_DIR_SRC ${CMAKE_SOURCE_DIR}/res)

target_include_directories(myapp PUBLIC
    "${RES_XYQ_DIR_SRC}"
)

# Create a variable that you can access directly during runtime as a preprocessor definition.
target_compile_definitions( myapp
    PRIVATE RES_PATH_XYQ="${RES_XYQ_DIR_SRC}"
)

# ======================================================================= <<<<<


# ======================================================================= >>>>>
#   Sub-Directories Inclusion
# ======================================================================= >>>>>


# ----------------------------------------------------------------- >>>>>
#   Utility (Personal) Inclusion
# ----------------------------------------------------------------- >>>>>

# Add a sub directory.
add_subdirectory( utils_XYQ )
target_link_libraries( myapp PUBLIC utils_XYQ )
target_include_directories(myapp PUBLIC
    "${PROJECT_SOURCE_DIR}/utils_XYQ"
)

# ----------------------------------------------------------------- <<<<<

# ======================================================================= <<<<<




