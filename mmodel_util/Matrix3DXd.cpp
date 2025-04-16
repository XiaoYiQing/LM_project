#include "Matrix3DXd.h"



// ====================================================================== >>>>>
//      Static Functions
// ====================================================================== >>>>>

bool Matrix3DXd::consist_check( const Matrix3DXd& tarMat ){
    return consist_check( tarMat.Mat3D );
}

bool Matrix3DXd::consist_check( const vector< Eigen::MatrixXd >& tarVec ){

    unsigned int lvl_cnt = tarVec.size();
    // The empty matrix case returns true always.
    if( lvl_cnt == 0 ){
        return true;
    }

    // Use the first matrix 
    unsigned int row_cnt_ref = tarVec.at(0).rows();
    unsigned int col_cnt_ref = tarVec.at(0).cols();

    // Matrix rows and cols consistency check.
    for( unsigned int z = 0; z < lvl_cnt; z++ ){
        if( tarVec.at(z).rows() != row_cnt_ref || tarVec.at(z).cols() != col_cnt_ref ){
            return false;
        }
    }

    // Reaching this point means no problem found.
    return true;

}


bool Matrix3DXd::same_size( const Eigen::MatrixXd& matA, const Eigen::MatrixXd& matB ){
    return ( matA.rows() == matB.rows() && matA.cols() == matB.cols() );
}


bool Matrix3DXd::null_ref_check( const Matrix3DXd& tarMat ){

    return null_ref_check( tarMat.Mat3D );

}

bool Matrix3DXd::null_ref_check( const vector< Eigen::MatrixXd >& tarVec ){

    return ( tarVec.at(0).rows() == 0 || tarVec.at(0).cols() == 0 );

}

// ====================================================================== <<<<<



// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

Matrix3DXd::Matrix3DXd(){
}

Matrix3DXd::Matrix3DXd( vector< Eigen::MatrixXd > Mat3D ){

    if( Mat3D.size() == 0 ){
        throw std::out_of_range( "Cannot initialize Matrix3DXd with an empty matrix vector." );
    }
    
    if( null_ref_check( Mat3D ) ){
        throw std::out_of_range( "Cannot initialize Matrix3DXd with an empty matrix vector." );
    }

    if( Matrix3DXd::consist_check( Mat3D ) ){
        this->Mat3D = Mat3D;
    }else{
        throw std::invalid_argument( "Matrix3DXd can only be initialized with a vector of Eigen::MatrixXd matrices with consistent dimensions." );
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
    size_vec[1] = this->Mat3D.at(0).cols();
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


bool Matrix3DXd::mat3DValidInCheck( const Eigen::MatrixXd& val ){

    // An empty 2D matrix is automatically invalidated.
    if( val.rows() == 0 || val.cols() == 0 ){
        return false;
    }

    // The empty matrix vector case returns true (Any non-trivial 2D mat can be added 
    // to an empty vec).
    if( this->Mat3D.size() == 0 ){
        return true;
    }

    // If neither input or mat vector are empty, just check for consistent 2D matrix dimensions.
    if( !same_size( val, this->Mat3D.at(0) ) ){
        return false;
    }

    // Reaching this point, no problem was found.
    return true;

}


void Matrix3DXd::push_back( const Eigen::MatrixXd& in2DMat ){

    if( this->mat3DValidInCheck( in2DMat ) ){
        this->Mat3D.push_back( in2DMat );
    }else{
        throw std::invalid_argument( "Matrix3DXd can only add 2D matrices with the same dimensions as the 2D matrices already present." );
    }

}

void Matrix3DXd::insert( unsigned int idx, const Eigen::MatrixXd& in2DMat ){

    if( idx > this->levels() ){
        throw std::out_of_range( "Insert index invalid: index is larger than the size of the vector." );
    }

    if( this->mat3DValidInCheck( in2DMat ) ){
        vector< Eigen::MatrixXd >::iterator iter = this->Mat3D.begin();
        iter = this->Mat3D.insert( iter + idx, in2DMat );
    }else{
        throw std::invalid_argument( "Matrix3DXd can only add 2D matrices with the same dimensions as the 2D matrices already present." );
    }
    
}

// ====================================================================== <<<<<
