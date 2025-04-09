cmake_minimum_required(VERSION 3.5)
project(MyFuncs)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


option(BUILD_SHARED_LIBS "Build shared libraries" OFF)


# Create our main executable.
add_executable( myapp src/main.cpp )




# ======================================================================= >>>>>
#   Sub-Directories Inclusion
# ======================================================================= >>>>>

# ======================================================================= <<<<<




