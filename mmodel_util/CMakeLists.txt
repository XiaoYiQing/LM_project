
# Create a list of variables, where the entries are the c++ and header files within the
# current directory.
set( mmodel_util_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/fData.h
	${CMAKE_CURRENT_SOURCE_DIR}/fData.cpp
)

add_library(mmodel_util ${mmodel_util_SOURCES})

target_link_libraries( mmodel_util PRIVATE Eigen3::Eigen )
