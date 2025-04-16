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

    
    

    Matrix3DXd_test_1(5);

    return 0; 

}



void Matrix3DXd_test_1( int case_idx ){

    // Initialize test case index.
    int case_cnt = 0;


    // 0- Empty vector initialization exception case.
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
    // 1- Multiplication operator test.
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
    // 2- Inplace multiplication operator test.
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


    case_cnt++;
    // 3- Invalid input vector initialization check.
    if( case_cnt == case_idx ){

        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }
        // Add another matrix to the vector, but having different dimensions 
        // than previous matrices.
        Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 3 );
        tmpMat << 9,9,9,9,9,9;
        tmp_mat_vec.push_back( tmpMat );

        try{
            Matrix3DXd my_3D_mat = Matrix3DXd( tmp_mat_vec );
        }catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        }

    }


    case_cnt++;
    // 4- push_back function check.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }
        
        // Use the vector of 2D matrices to initialize a 3D matrix.
        Matrix3DXd my_3D_mat;
        try{
            my_3D_mat = Matrix3DXd( tmp_mat_vec );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        }

        // Add another matrix to the vector at its end.
        Eigen::MatrixXd addMat = Eigen::MatrixXd( 2, 2 );
        addMat << 91, 92, 93, 94;
        try{
            my_3D_mat.push_back( addMat );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        }
        cout << my_3D_mat.at( my_3D_mat.levels() - 1 ) << endl;

    }



    case_cnt++;
    // 5- insert function check.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }

        // Use the vector of 2D matrices to initialize a 3D matrix.
        Matrix3DXd my_3D_mat;
        try{
            my_3D_mat = Matrix3DXd( tmp_mat_vec );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        }
        
        // Insert another matrix to the vector at the target valid index.
        Eigen::MatrixXd addMat = Eigen::MatrixXd( 2, 2 );
        addMat << 91, 92, 93, 94;
        unsigned int tarIdx = 2;
        try{
            my_3D_mat.insert( tarIdx, addMat );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        } catch( const std::out_of_range& e ){
            cerr << e.what() << endl;
            return;
        }
        cout << my_3D_mat.at( tarIdx ) << endl;
        cout << my_3D_mat.at( tarIdx - 1 ) << endl;


        // Insert another matrix to the vector at the target invalid index.
        Eigen::MatrixXd addMat2 = Eigen::MatrixXd( 2, 2 );
        addMat2 << 51, 52, 53, 54;
        tarIdx = 100;
        try{
            my_3D_mat.insert( tarIdx, addMat2 );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        } catch( const std::out_of_range& e ){
            cerr << e.what() << endl;
            return;
        }
        cout << my_3D_mat.at( tarIdx ) << endl;
        cout << my_3D_mat.at( tarIdx - 1 ) << endl;

    }


    return;

}


