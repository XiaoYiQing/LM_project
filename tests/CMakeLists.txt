
# Create a list of variables, where the entries are the c++ and header files within the
# current directory.
set( tests_SOURCES 
	${CMAKE_CURRENT_SOURCE_DIR}/tests_Eigen.h
	${CMAKE_CURRENT_SOURCE_DIR}/tests_Eigen.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tests_Matrix3DXd.h
	${CMAKE_CURRENT_SOURCE_DIR}/tests_Matrix3DXd.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tests_fData.h
	${CMAKE_CURRENT_SOURCE_DIR}/tests_fData.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tests_numUtils.h
	${CMAKE_CURRENT_SOURCE_DIR}/tests_numUtils.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tests_LM_eng.h
	${CMAKE_CURRENT_SOURCE_DIR}/tests_LM_eng.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tests_LM_eng_2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tests_LTI_descSyst.h
	${CMAKE_CURRENT_SOURCE_DIR}/tests_LTI_descSyst.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tests_orchestrators.h
	${CMAKE_CURRENT_SOURCE_DIR}/tests_orchestrators.cpp
)

add_library(tests ${tests_SOURCES})

target_link_libraries( tests PRIVATE Eigen3::Eigen )

target_link_libraries( tests PRIVATE mmodel_util )
target_include_directories(tests PUBLIC
    "${PROJECT_SOURCE_DIR}/mmodel_util"
)

target_link_libraries( tests PRIVATE utils_XYQ )
target_include_directories(tests PUBLIC
    "${PROJECT_SOURCE_DIR}/utils_XYQ"
)
