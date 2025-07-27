#ifndef MATRIX3DXCD_H
#define MATRIX3DXCD_H


#include <Eigen/Dense>
#include <iostream>
#include <magic_enum.hpp>
#include <string>   
#include <vector> 


using namespace std;



/*
This is my take for the 3D matrix based on the MatrixXcd matrix from the 
Eigen library.
The 3D matrix is nothing more than a vector of 2D MatrixXcd matrices.
The class' purpose is to ease manipulation of the 3D matrix as well as simple
matrix operations.
This implementation doesn't operate as a typical 3D matrix in that its 3rd dimension
is variable.
*/
class Matrix3DXcd{    


public:

// ====================================================================== >>>>>
//      Static Functions
// ====================================================================== >>>>>

    static const double DEF_NUM_THRESH;

    /*
    Function for checking the consistency of the 3D matrix vector.
    For instance, all 2D matrices in the vector must have the same size.
    */
    static bool consist_check( const Matrix3DXcd& );
    /*
    Function for checking the consistency of the 3D matrix vector.
    For instance, all 2D matrices in the vector must have the same size.
    */
    static bool consist_check( const vector< Eigen::MatrixXcd >& );

    /*
    Function for checking if 2 2D matrices are of the same size.
    */
    static bool same_size( const Eigen::MatrixXcd& matA, const Eigen::MatrixXcd& matB );

    /*
    Function for checking if the reference matrix of the matrix vector has
    0 rows or 0 columns.
    If true, the reference matrix has 0 rows or columns.
    */
    static bool null_ref_check( const Matrix3DXcd& );
    /*
    Function for checking if the reference matrix of the matrix vector has
    0 rows or 0 columns.
    If true, the reference matrix has 0 rows or columns.
    */
    static bool null_ref_check( const vector< Eigen::MatrixXcd >& );


// ====================================================================== <<<<<



protected:

    
// ====================================================================== >>>>>
//      Member Variables
// ====================================================================== >>>>>

    /* 
    The numerical threshold used by this object to determine if a numerical value
    is to be considered trivial.
    As such, this value is also used to evaluate if two numerical values are equal.
    */
    double num_thresh;

    /*
    The 3D matrix, which is just a vector of 2D matrices.
    */
    vector< Eigen::MatrixXd > Mat3D;

// ====================================================================== <<<<<

};



#endif  // MATRIX3DXCD_H