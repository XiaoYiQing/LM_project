#include "tests_Matrix3DXcd.h"



void tests::Matrix3DXcd_serialize_test(){

    // Create a vector of 2D matrices.
    unsigned int level_cnt = 10;
    vector< Eigen::MatrixXcd > tmp_mat_vec;
    for( unsigned int z = 0; z < level_cnt; z++ ){
        Eigen::MatrixXcd tmpMat = Eigen::MatrixXcd( 2, 2 );
        tmpMat << 0+z*10, 1+z*10, 2+z*10, 3+z*10;
        tmp_mat_vec.push_back( tmpMat );
    }
    // Initialize a 3D matrix using the vector.
    Matrix3DXcd myMat3D = Matrix3DXcd( tmp_mat_vec );

    string outFileFullFileName = SRC_PATH_XYQ_str + "/data_output/Matrix3DXcd.bin";

    myMat3D.serialize( outFileFullFileName );

    Matrix3DXcd myMat3D_2;

    myMat3D_2.deserialize( outFileFullFileName );

    bool test_bool = true;
    for( unsigned int z = 0; z < level_cnt; z++ ){
        test_bool = test_bool && ( myMat3D.at(z) == myMat3D_2.at(z) );
    }

    if( test_bool ){
        cout << "Matrix3DXcd serialize/deserialize test: Passed!" << endl;
    }else{
        cout << "Matrix3DXcd serialize/deserialize test: failed!" << endl;
    }

}