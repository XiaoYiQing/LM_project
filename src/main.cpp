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


using namespace std;

// Global string variable holding the full path of the extra resource directory.
string RES_PATH_XYQ_str = string( RES_PATH_XYQ );


void Matrix3DXd_test_1( int case_idx );


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

    
    

    Matrix3DXd_test_1(2);

    return 0; 

}



void Matrix3DXd_test_1( int case_idx ){

    // Initialize test case index.
    int case_cnt = 0;


    // Empty vector initialization exception case.
    if( case_cnt == case_idx ){

        vector< Eigen::MatrixXd > tmp_mat_vec;

        try{
            Matrix3DXd my_3D_mat = Matrix3DXd( tmp_mat_vec );
        }catch (const std::out_of_range& e){
            cerr << e.what() << endl;
            return;
        }

    }

    case_cnt++;
    // Multiplication operator test.
    if( case_cnt == case_idx ){

        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }

        Matrix3DXd my3DMat = Matrix3DXd( tmp_mat_vec );
        cout << my3DMat.at(2) << endl;
        Matrix3DXd my3DMat2 = my3DMat*3;
        cout << my3DMat2.at(2) << endl;
                
    }
    
    case_cnt++;
    // Inplace multiplication operator test.
    if( case_cnt == case_idx ){

        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }

        Matrix3DXd my3DMat = Matrix3DXd( tmp_mat_vec );
        cout << my3DMat.at(2) << endl;
        my3DMat *= 3;
        cout << my3DMat.at(2) << endl;
                
    }


    return;

}


