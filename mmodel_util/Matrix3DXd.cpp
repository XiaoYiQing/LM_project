#include "Matrix3DXd.h"



// ====================================================================== >>>>>
//      Static Functions
// ====================================================================== >>>>>

bool Matrix3DXd::consist_check( const vector< Eigen::MatrixXd >& tarMat ){

    unsigned int lvl_cnt = tarMat.size();
    // The empty matrix case returns true always.
    if( lvl_cnt == 0 ){
        return true;
    }

    // Use the first matrix 
    unsigned int row_cnt_ref = tarMat.at(0).rows();
    unsigned int col_cnt_ref = tarMat.at(0).cols();

    // Matrix rows and cols consistency check.
    for( unsigned int z = 0; z < lvl_cnt; z++ ){
        if( tarMat.at(z).rows() != row_cnt_ref || tarMat.at(z).cols() != col_cnt_ref ){
            return false;
        }
    }

    // Reaching this point means no problem found.
    return true;

}

// ====================================================================== <<<<<



// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

Matrix3DXd::Matrix3DXd( unsigned int row_cnt, unsigned int col_cnt ){
    this->row_cnt = row_cnt;
    this->col_cnt = col_cnt;
}

Matrix3DXd::Matrix3DXd( vector< Eigen::MatrixXd > Mat3D ){

    if( Mat3D.size() == 0 ){
        throw std::out_of_range( "Cannot initialize Matrix3DXd with an empty matrix vector." );
    }
    
    if( Matrix3DXd::consist_check( Mat3D ) ){
        this->Mat3D = Mat3D;
        this->row_cnt = Mat3D.at(0).rows();
        this->col_cnt = Mat3D.at(0).rows();
    }else{
        throw std::invalid_argument( "Matrix3DXd can only be initialized with a vector of Eigen::MatrixXd matrices sharing the same dimensions." );
    }

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
    for( unsigned int i = 0; i < lvl_cnt; i++ ){
        Mat3D_res.Mat3D[i] *= scalar;
    }
    
    return Mat3D_res;
}

Matrix3DXd& Matrix3DXd::operator*=(const double scalar){

    // Obtain 3D matrix level count.
    unsigned int lvl_cnt = this->Mat3D.size();
    // Apply individual 3D matrix scalar multiplication.
    for( unsigned int i = 0; i < lvl_cnt; i++ ){
        this->Mat3D[i] *= scalar;
    }

    return *this;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>


vector<unsigned int> Matrix3DXd::size() const{

    vector<unsigned int> size_vec = {0,0,0};
    if( this->Mat3D.size() == 0 ){
        return size_vec;
    }

    size_vec[0] = this->Mat3D.at(0).rows();
    size_vec[1] = this->Mat3D.at(1).cols();
    size_vec[2] = this->Mat3D.size();

    return size_vec;

}

unsigned int Matrix3DXd::rows() const{
    if( Mat3D.size() == 0 ){
        return 0;
    }else{
        return Mat3D.at(0).rows();
    }
}
unsigned int Matrix3DXd::cols() const{
    if( Mat3D.size() == 0 ){
        return 0;
    }else{
        return Mat3D.at(0).cols();
    }
}
unsigned int Matrix3DXd::levels() const{
    return Mat3D.size();
}

bool Matrix3DXd::isEmpty(){
    return this->Mat3D.size() == 0;
}

Eigen::MatrixXd Matrix3DXd::at( unsigned int z ) const{
    return this->Mat3D.vector::at(z);
}


void Matrix3DXd::push_back( const Eigen::MatrixXd& val ){

    if( this->isEmpty() ){
        this->Mat3D.push_back( val );
        return;
    }

    if( val.rows() == this->Mat3D.at(0).rows() && val.cols() == this->Mat3D.at(0).cols() ){
        

    }else{

        
    }

}

// ====================================================================== <<<<<
