#include "Matrix3DXd.h"




// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

Matrix3DXd::Matrix3DXd(){

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
