#include "Matrix3DXd.h"




// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

Matrix3DXd::Matrix3DXd(){

}

Matrix3DXd::Matrix3DXd( vector< Eigen::MatrixXd > Mat3D ){
    this->Mat3D = Mat3D;
}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Operators
// ====================================================================== >>>>>

// Multiplication operator overload.
Matrix3DXd Matrix3DXd::operator*(const double scalar) const {

    // Obtain 3D matrix level count.
    unsigned int lvl_cnt = Mat3D.size();
    // New 3D matrix definition.
    Matrix3DXd Mat3D_res = Matrix3DXd( this->Mat3D );
    // Apply individual 3D matrix scalar multiplication.
    for( int i = 0; i < lvl_cnt; i++ ){
        Mat3D_res.Mat3D[i] *= scalar;
    }
    
    return Mat3D_res;
}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>


vector<unsigned int> Matrix3DXd::size(){

    vector<unsigned int> size_vec = {0,0,0};
    if( Mat3D.size() == 0 ){
        return size_vec;
    }

    size_vec[0] = Mat3D.at(0).rows();
    size_vec[1] = Mat3D.at(1).cols();
    size_vec[2] = Mat3D.size();

}

unsigned int Matrix3DXd::rows(){
    if( Mat3D.size() == 0 ){
        return 0;
    }else{
        return Mat3D.at(0).rows();
    }
}
unsigned int Matrix3DXd::cols(){
    if( Mat3D.size() == 0 ){
        return 0;
    }else{
        return Mat3D.at(0).cols();
    }
}
unsigned int Matrix3DXd::levels(){
    return Mat3D.size();
}



// ====================================================================== <<<<<
