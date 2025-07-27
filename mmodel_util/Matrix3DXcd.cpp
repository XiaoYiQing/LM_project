#include "Matrix3DXcd.h"



// ====================================================================== >>>>>
//      Static Functions
// ====================================================================== >>>>>

/*
Default numerical threshold for determining if a number is trivial.
*/
const double Matrix3DXcd::DEF_NUM_THRESH = 1e-12;

bool Matrix3DXcd::consist_check( const Matrix3DXcd& tarMat ){
    return consist_check( tarMat.Mat3D );
}

bool Matrix3DXcd::consist_check( const vector< Eigen::MatrixXcd >& tarVec ){

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


bool Matrix3DXcd::same_size( const Eigen::MatrixXcd& matA, const Eigen::MatrixXcd& matB ){
    return ( matA.rows() == matB.rows() && matA.cols() == matB.cols() );
}


bool Matrix3DXcd::null_ref_check( const Matrix3DXcd& tarMat ){

    return null_ref_check( tarMat.Mat3D );

}

bool Matrix3DXcd::null_ref_check( const vector< Eigen::MatrixXcd >& tarVec ){

    return ( tarVec.at(0).rows() == 0 || tarVec.at(0).cols() == 0 );

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

Matrix3DXcd::Matrix3DXcd(){
    this->num_thresh = Matrix3DXcd::DEF_NUM_THRESH;
}

Matrix3DXcd::Matrix3DXcd( unsigned int row_idx, unsigned int col_idx, unsigned int lvl_idx ){
    this->num_thresh = Matrix3DXcd::DEF_NUM_THRESH;

    this->reInit( row_idx, col_idx, lvl_idx );
}

Matrix3DXcd::Matrix3DXcd( vector< Eigen::MatrixXcd > Mat3D ){

    this->num_thresh = Matrix3DXcd::DEF_NUM_THRESH;

    if( Mat3D.size() == 0 ){
        throw std::out_of_range( "Cannot initialize Matrix3DXcd with an empty matrix vector." );
    }
    
    if( null_ref_check( Mat3D ) ){
        throw std::out_of_range( "Cannot initialize Matrix3DXcd with an empty matrix vector." );
    }

    if( Matrix3DXcd::consist_check( Mat3D ) ){
        this->Mat3D = Mat3D;
    }else{
        throw std::invalid_argument( "Matrix3DXcd can only be initialized with a vector of Eigen::MatrixXcd matrices with consistent dimensions." );
    }


}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Operations
// ====================================================================== >>>>>

void Matrix3DXcd::elem_pow( double exp_val ){
    for( unsigned int z = 0; z < Mat3D.size(); z++ ){
        this->Mat3D.at(z) = this->Mat3D.at(z).array().pow( exp_val );
    }
}

void Matrix3DXcd::elem_raise_pow( double base_val ){
    
    double log_base_val = log( base_val );

    *this*=log_base_val;

    unsigned int z = 0;
    for( unsigned int z = 0; z < Mat3D.size(); z++ ){
        this->Mat3D.at(z) = this->Mat3D.at(z).array().exp();
    }

}

void Matrix3DXcd::elem_log10(){
    for( unsigned int z = 0; z < Mat3D.size(); z++ ){
        this->Mat3D.at(z) = this->Mat3D.at(z).array().log10();
    }
}

void Matrix3DXcd::elem_cos(){
    for( unsigned int z = 0; z < Mat3D.size(); z++ ){
        this->Mat3D.at(z) = this->Mat3D.at(z).array().cos();
    }
}

void Matrix3DXcd::elem_sin(){
    for( unsigned int z = 0; z < Mat3D.size(); z++ ){
        this->Mat3D.at(z) = this->Mat3D.at(z).array().sin();
    }
}

void Matrix3DXcd::elem_atan(){
    for( unsigned int z = 0; z < Mat3D.size(); z++ ){
        this->Mat3D.at(z) = this->Mat3D.at(z).array().atan();
    }
}

Matrix3DXcd Matrix3DXcd::elem_div_spec( const Matrix3DXcd& tarMat ) const{

    unsigned int currLevels = this->levels();

    if( currLevels != tarMat.levels() ){
        throw std::invalid_argument( "Element-wise multipying matrix has mismatched dimensions." );
    }
    if( !Matrix3DXcd::consist_check( tarMat ) ){
        throw std::invalid_argument( "Element-wise multipying matrix has inconsistent dimensions." );
    }
    if( !Matrix3DXcd::same_size( tarMat.at(0), this->Mat3D.at(0) ) ){
        throw std::invalid_argument( "Element-wise multipying matrix has mismatched dimensions." );
    }

    unsigned int row_cnt = this->rows();
    unsigned int col_cnt = this->cols();

    Matrix3DXcd resMat;
    resMat.reInit( row_cnt, col_cnt, currLevels );

    for( unsigned int z = 0; z < currLevels; z++ ){
        // A more careful evaluation 
        for( unsigned int i = 0; i < row_cnt; i++ ){
            for( unsigned int j = 0; j < col_cnt; j++ ){

                if( abs( Mat3D.at(z)(i,j) ) < this->num_thresh ){
                    resMat.Mat3D.at(z)(i,j) = 0;
                }
                else if( abs( tarMat.at(z)(i,j) ) >= this->num_thresh ){
                    resMat.Mat3D.at(z)(i,j) = Mat3D.at(z)(i,j)/tarMat.at(z)(i,j);
                }
                else{
                    throw std::runtime_error( "Division by zero error." );
                }

            }
        }
    }

    return resMat;
}

// ====================================================================== <<<<<



// ====================================================================== >>>>>
//      Specialized Operations
// ====================================================================== >>>>>

Matrix3DXcd Matrix3DXcd::elem_phase_comp( const Matrix3DXcd& tarMat ){

    unsigned int reLevels = tarMat.levels();
    unsigned int row_cnt = tarMat.rows();
    unsigned int col_cnt = tarMat.cols();

    Matrix3DXcd resMat;
    resMat.reInit( row_cnt, col_cnt, reLevels );

    for( unsigned int z = 0; z < reLevels; z++ ){
    for( unsigned int i = 0; i < row_cnt; i++ ){
    for( unsigned int j = 0; j < col_cnt; j++ ){
        resMat.Mat3D.at(z)(i,j) = atan( tarMat.Mat3D.at(z)(i,j) );
    }    }    }

    return resMat;

}


double Matrix3DXcd::RMS_total_comp( const Matrix3DXcd& tarMat ){

    unsigned int reLevels = tarMat.levels();
    unsigned int row_cnt = tarMat.rows();
    unsigned int col_cnt = tarMat.cols();

    double resVal = 0;
    double currRe = 0, currIm = 0;

    // Compute the sum of all values squared.
    for( unsigned int z = 0; z < reLevels; z++ ){
    for( unsigned int i = 0; i < row_cnt; i++ ){
    for( unsigned int j = 0; j < col_cnt; j++ ){
        resVal += pow( tarMat.Mat3D.at(z)(i,j).real(), 2 ) + 
            pow( tarMat.Mat3D.at(z)(i,j).imag(), 2 );
    }    }    }
    // Compute the mean square value.
    resVal /= ( reLevels * row_cnt * col_cnt );
    // Take the square root of the final mean square value.
    resVal = pow( resVal, 0.5 );

    return resVal;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Support Functions
// ====================================================================== >>>>>

bool Matrix3DXcd::mat3DValidInCheck( const Eigen::MatrixXcd& val ){

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

void Matrix3DXcd::reInit( unsigned int row_cnt, unsigned int col_cnt, unsigned int lvl_cnt ){

    if( row_cnt == 0 || col_cnt == 0 || lvl_cnt == 0 ){
        throw std::domain_error( "2D matrix vector cannot reserve if its 2D matrices have 0 row or column." );
    }

    this->Mat3D.clear();

    Eigen::MatrixXcd firstMat = Eigen::MatrixXcd( row_cnt, col_cnt );
    firstMat.setZero();
    this->push_back( firstMat );
    this->reserve( lvl_cnt );

}

void Matrix3DXcd::clear(){
    this->Mat3D.clear();
}

void Matrix3DXcd::reserve( unsigned int res_size ){

    size_t currSize = this->levels();
    unsigned int row_cnt = this->rows();
    unsigned int col_cnt = this->cols();

    if( row_cnt == 0 || col_cnt == 0 ){
        throw std::invalid_argument( "2D matrix vector cannot reserve if its 2D matrices have 0 row or column." );
    }

    this->Mat3D.reserve( res_size );
    for( size_t i = currSize; i < res_size; i++ ) {
        Mat3D.emplace_back( Eigen::MatrixXcd( row_cnt, col_cnt ) );
    }

}


void Matrix3DXcd::resize( unsigned int tarSize ){
    this->Mat3D.resize( tarSize );
}

void Matrix3DXcd::shrink_to_fit(){
    this->Mat3D.shrink_to_fit();
}

Matrix3DXcd Matrix3DXcd::segment( size_t startIndex, size_t len ) const{

    // Obtain the continuous subset of 2D matrices using iterators.
    vector<Eigen::MatrixXcd> subset(this->Mat3D.begin() + startIndex, 
        this->Mat3D.begin() + startIndex + len);

    Matrix3DXcd segMat( subset );

    return segMat;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>


vector<unsigned int> Matrix3DXcd::size() const{

    vector<unsigned int> size_vec = {0,0,0};
    if( this->Mat3D.size() == 0 ){
        return size_vec;
    }

    size_vec[0] = this->Mat3D.at(0).rows();
    size_vec[1] = this->Mat3D.at(0).cols();
    size_vec[2] = this->Mat3D.size();

    return size_vec;

}

unsigned int Matrix3DXcd::rows() const{
    if( Mat3D.size() == 0 ){
        return 0;
    }else{
        return Mat3D.at(0).rows();
    }
}
unsigned int Matrix3DXcd::cols() const{
    if( Mat3D.size() == 0 ){
        return 0;
    }else{
        return Mat3D.at(0).cols();
    }
}
unsigned int Matrix3DXcd::levels() const{
    return Mat3D.size();
}

bool Matrix3DXcd::isEmpty() const{
    return this->Mat3D.size() == 0;
}

Eigen::MatrixXcd Matrix3DXcd::at( unsigned int z ) const{
    return this->Mat3D.vector::at(z);
}

Matrix3DXcd Matrix3DXcd::at( vector< unsigned int > idxVec ){
    
    vector<int> idxVec_tmp;
    idxVec_tmp.reserve( idxVec.size() );
    for( unsigned int z : idxVec ){
        idxVec_tmp.push_back( z );
    }

    Eigen::VectorXi idxVecAlt = Eigen::Map<Eigen::VectorXi>(idxVec_tmp.data(), idxVec_tmp.size());

    // Check for out of range violation.
    if( (unsigned int) idxVecAlt.maxCoeff() > this->levels() ){
        throw std::out_of_range( "Index array contain at least one out of range index." );
    }

    // Initialize sub-matrix vector.
    vector< Eigen::MatrixXcd > subMatVec;
    // Reserve memory for sub-vector.
    subMatVec.reserve( idxVecAlt.size() );
    // Fill the sub-vector with the targeted entries.
    for( unsigned int z : idxVecAlt ){
        subMatVec.emplace_back( this->Mat3D.at(z) );
    }

    // Return the vector under the form of a 3D matrix.
    return Matrix3DXcd( subMatVec );

}

void Matrix3DXcd::set( unsigned int tarIdx, const Eigen::MatrixXcd& tarMat ){

    if( tarIdx >= this->levels() ){
        throw std::out_of_range( 
            "Specified vector index outside range of available vector entries." );
    }

    if( !same_size( tarMat, this->Mat3D.at(0) ) ){
        throw std::invalid_argument( 
            "Input 2D matrix must have the same dimensions as the 2D matrices already present." );
    }

    // Assign the new 2D matrix values.
    this->Mat3D.at( tarIdx ) = tarMat;

}


void Matrix3DXcd::set( unsigned int row_idx, unsigned int col_idx, unsigned int lvl_idx, double val ){

    if( lvl_idx >= this->levels() ){
        throw std::out_of_range( 
            "Specified vector index outside range of available vector entries." );
    }

    if( row_idx >= this->rows() || col_idx >= this->cols() ){
        throw std::out_of_range( 
            "2D matrix coordinate is out of range." );
    }

    // Assign the new 2D matrix values.
    Mat3D.at( lvl_idx )( row_idx, col_idx ) = val;

}




void Matrix3DXcd::push_back( const Eigen::MatrixXcd& in2DMat ){

    if( this->mat3DValidInCheck( in2DMat ) ){
        this->Mat3D.push_back( in2DMat );
    }else{
        throw std::invalid_argument( "Matrix3DXcd can only add 2D matrices with the same dimensions as the 2D matrices already present." );
    }

}

void Matrix3DXcd::push_back( const vector< Eigen::MatrixXcd >& in2DMatVec ){

    unsigned int inSize = in2DMatVec.size();

    // Empty vector check.
    if( inSize == 0 ){
        return;
    }

    // Check dimension consistency of matrices within the vector.
    if( !Matrix3DXcd::consist_check( in2DMatVec ) ){
        throw std::invalid_argument( "Matrix3DXcd::push_back: \n\tInput vector of matrices must have consistent matrix dimensions." );
    }

    // Check for empty 2D matrix special case.
    if( this->mat3DValidInCheck( in2DMatVec.at(0) ) ){
        // Reserve the required memory.
        this->Mat3D.reserve( this->Mat3D.size() + inSize );
        // Put all matrices at the end of the vector in the given order.
        for( unsigned int z = 0; z < inSize; z++ ){
            Mat3D.push_back( in2DMatVec.at(z) );
        }

    }else{
        throw std::invalid_argument( "Matrix3DXcd::push_back: \n\tInput matrices are invalid." );
    }


}

void Matrix3DXcd::insert( unsigned int idx, const Eigen::MatrixXcd& in2DMat ){

    // Index validity check.
    if( idx > this->levels() ){
        throw std::out_of_range( "Insert index invalid: index is larger than the size of the vector." );
    }

    // Attempt at insert.
    if( this->mat3DValidInCheck( in2DMat ) ){
        vector< Eigen::MatrixXcd >::iterator iter = this->Mat3D.begin();
        iter = this->Mat3D.insert( iter + idx, in2DMat );
    }else{
        throw std::invalid_argument( "Matrix3DXcd can only add 2D matrices with the same dimensions as the 2D matrices already present." );
    }
    
}

// ====================================================================== <<<<<
