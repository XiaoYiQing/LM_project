#ifndef MATRIX3DXD_H
#define MATRIX3DXD_H


#include <Eigen/Dense>
#include <magic_enum.hpp>
#include <string>   
#include <vector> 



using namespace std;



/*
This is my take for the 3D matrix based on the MatrixXd matrix from the 
Eigen library.
The 3D matrix is nothing more than a vector of 2D MatrixXd matrices.
The class' purpose is to ease manipulation of the 3D matrix as well as simple
matrix operations.
This implementation doesn't operate as a typical 3D matrix in that its 3rd dimension
is variable.
*/
class Matrix3DXd{    


public:

// ====================================================================== >>>>>
//      Static Functions
// ====================================================================== >>>>>

    /*
    Function for checking the consistency of the 3D matrix vector.
    For instance, all 2D matrices in the vector must have the same size.
    */
    static bool consist_check( const Matrix3DXd& );
    /*
    Function for checking the consistency of the 3D matrix vector.
    For instance, all 2D matrices in the vector must have the same size.
    */
    static bool consist_check( const vector< Eigen::MatrixXd >& );

    /*
    Function for checking if 2 2D matrices are of the same size.
    */
    static bool same_size( const Eigen::MatrixXd& matA, const Eigen::MatrixXd& matB );

    /*
    Function for checking if the reference matrix of the matrix vector has
    0 rows or 0 columns.
    If true, the reference matrix has 0 rows or columns.
    */
    static bool null_ref_check( const Matrix3DXd& );
    /*
    Function for checking if the reference matrix of the matrix vector has
    0 rows or 0 columns.
    If true, the reference matrix has 0 rows or columns.
    */
    static bool null_ref_check( const vector< Eigen::MatrixXd >& );


    

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

    Matrix3DXd();

    Matrix3DXd( vector< Eigen::MatrixXd > Mat3D );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Operators
// ====================================================================== >>>>>

    // Overloading the multiplication operator for scalar multiplication.
    Matrix3DXd operator*(const double scalar) const;

    // Overloading the compound multiplication operator for scalar multiplication.
    Matrix3DXd& operator*=(const double scalar);

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Operations
// ====================================================================== >>>>>

    /*
    Perform element-wise power to the target value on the entries of the matrix.
    */
    void elem_pow( double exp_val );

    /*
    Perform element-wise power with the input base value raised to the power value given
    by the entries of the matrix.
    */
    void elem_raise_pow( double base_val );

    /*
    Perform element-wise base 10 logarithmic operation.
    */
    void elem_log10(  );

// ====================================================================== <<<<<



// ====================================================================== >>>>>
//      Support Functions
// ====================================================================== >>>>>

    /* 
    Check if the target 2D matrix can be added to the 3D matrix.
    */
    bool mat3DValidInCheck( const Eigen::MatrixXd& val );

    /*
    Re-initialize the 3D matrix with an initiate set of zero matrices of specified sizes.
    This function effectively serves as the "reserve" function specifically when the 3D matrix
    is empty (or emptied) and requires a new start with specificed number of rows and columns.
    WARNING: This function erases all current 3D matrix entries.
    */
    void reInit( unsigned int row_idx, unsigned int col_idx, unsigned int lvl_idx );

    /*
    Reserve a set number of entries worth of memory for additional 2D matrices.
    This reserve additionally places zero matrices at the reserved spaces within
    the 3D matrix as placeholders for direct index access.
    */
    void reserve( unsigned int );

    /*
    Perform the same operations as resize() of a vector variable.
    */
    void resize( unsigned int );

    /*
    Perform the same operations as shrink_to_fit() of a vector variable.
    */
    void shrink_to_fit();

// ====================================================================== <<<<<



// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>

    // Obtain a vector representing the dimensions of the 3D matrix (row #, col #, level # ).
    vector<unsigned int> size() const;

    unsigned int rows() const;
    unsigned int cols() const;
    unsigned int levels() const;

    // Determine if the matrix is empty. True if matrix vector is empty.
    bool isEmpty();

    /*
    Return the 2D matrix at the target index of the vector of 2D matrix.
    NOTE: this "at( unsigned int )" function does not perform the assignment 
    function because it doesn't return a reference, but a copy. This is to 
    prevent insertion of a 2D matrix having different dimensions than the rest.
    */
    Eigen::MatrixXd at( unsigned int ) const;

    /*
    Set the 2D matrix values at the target index to the given 2D matrix values.
    */
    void set( unsigned int, const Eigen::MatrixXd& );

    /*
    Set the value at the specific coordinate of the specific 2D matrix. 
    */
    void set( unsigned int row_idx, unsigned int col_idx, unsigned int lvl_idx, double val );

    /*
    Put a new 2D matrix entry at the end of the vector.
    */
    void push_back( const Eigen::MatrixXd& in2DMat );
    /*
    Put a vector of 2D matrix entries at the end of the vector.
    */
    void push_back( const vector< Eigen::MatrixXd >& in2DMatVec );

    /*
    Insert a new 2D matrix entry at the target index position.
    */
    void insert( unsigned int idx, const Eigen::MatrixXd& in2DMat );


// ====================================================================== <<<<<






protected:

    
// ====================================================================== >>>>>
//      Member Variables
// ====================================================================== >>>>>

    /*
    The 3D matrix, which is just a vector of 2D matrices.
    */
    vector< Eigen::MatrixXd > Mat3D;

// ====================================================================== <<<<<

private:




};





#endif  // MATRIX3DXD_H