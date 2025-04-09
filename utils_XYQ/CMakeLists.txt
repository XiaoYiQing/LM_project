
# Create a list of variables, where the entries are the c++ and header files within the
# current directory.
set( utils_XYQ_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/timeUtils.h
	${CMAKE_CURRENT_SOURCE_DIR}/timeUtils.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/numUtils.h
	${CMAKE_CURRENT_SOURCE_DIR}/numUtils.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/magic_enum.hpp 
)

add_library(utils_XYQ ${utils_XYQ_SOURCES})
