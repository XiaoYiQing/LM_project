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
#include "tests_Matrix3DXd.h"


using namespace std;

// Global string variable holding the full path of the extra resource directory.
string RES_PATH_XYQ_str = string( RES_PATH_XYQ );



int main() {

    int err_ret_val = 1;

    cout << RES_PATH_XYQ_str << endl;

    // Define our frequency data object.
    fData myF;

    // Test aspectes of the Eigen library.
    // int testCase = 3;
    // tests::eigen_test1( testCase );

    // Define the full file name.
    // string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
    // fData::read_sXp_file( myF, fullFileName );

    // cout << myF.get_f_cnt() << endl;
    // cout << myF.get_f_scale_str() << endl;
    // cout << myF.get_f_scale_num() << endl;
    // cout << myF.get_reData_at_f( 10 ) << endl;
    // cout << myF.get_imData_at_f( 10 ) << endl;


    // tests::Matrix3DXd_test_1(7);

    // Create the 3D matrix object using the above matrix vector.
    Matrix3DXd my_3D_mat = Matrix3DXd();

    Eigen::MatrixXd init_mat = Eigen::MatrixXd( 2, 2 );
    init_mat << 1, 2, 3, 4;
    my_3D_mat.push_back( init_mat );

    Eigen::MatrixXd next_mat = Eigen::MatrixXd( 2, 2 );
    next_mat << 11, 12, 13, 14;
    my_3D_mat.at(0) = next_mat;

    cout << my_3D_mat.at(0) << endl;

    return 0; 

}



