// Preprocessor variable created when running CMAKE process. This variable 
// indicates the location where non-compilabe or non c/c++ support files are located.
#ifndef RES_PATH_XYQ
#   define RES_PATH_XYQ ""
#else
#endif

#include <Eigen/Dense>
#include <iostream>
#include <filesystem>
#include <fstream>
#include "Matrix3DxD.h"
#include <regex>
#include <sstream>
#include <string>   
#include <vector>  



#include "fData.h"
#include "numUtils.h"
#include "tests_Eigen.h"
#include "tests_fData.h"
#include "tests_Matrix3DXd.h"
#include "tests_numUtils.h"


using namespace std;

// Global string variable holding the full path of the extra resource directory.
string RES_PATH_XYQ_str = string( RES_PATH_XYQ );



int main() {

    int err_ret_val = 1;

    

    // Test aspectes of the Eigen library.
    // int testCase = 3;
    // tests::eigen_test1( testCase );


    tests::numUtils_test_1(2);

    // tests::Matrix3DXd_test_1(9);
    // tests::Matrix3DXd_test_2_ops(3);
    // tests::Matrix3DXd_test_3_spec_ops(0);
    // tests::Matrix3DXd_test_4_supp(0);

    // tests::fData_test_1( 4 );
    // tests::fData_test_2( 1 );


    
    

    return 0; 

}



