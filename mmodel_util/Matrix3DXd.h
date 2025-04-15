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


    /*
    Function for checking the integrety of the 3D matrix vector.
    For instance, all 2D matrices must have the same size.
    */
    static bool consist_check( const vector< Eigen::MatrixXd >& );


// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

    Matrix3DXd( unsigned int row_cnt = 0, unsigned int col_cnt = 0 );

    Matrix3DXd( vector< Eigen::MatrixXd > Mat3D );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Operators
// ====================================================================== >>>>>

    // Overloading the multiplication operator
    Matrix3DXd operator*(const double scalar) const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>

    // Obtain a vector representing the dimensions of the 3D matrix (row #, col #, level # ).
    vector<unsigned int> size() const;

    unsigned int rows() const;
    unsigned int cols() const;
    unsigned int levels() const;

// ====================================================================== <<<<<






protected:

    
// ====================================================================== >>>>>
//      Member Variables
// ====================================================================== >>>>>
    
    // Number of rows.
    unsigned int row_cnt;
    // Number of columns.
    unsigned int col_cnt;

    /*
    The 3D matrix, which is just a vector of 2D matrices.
    */
    vector< Eigen::MatrixXd > Mat3D;

// ====================================================================== <<<<<

private:




};





#endif  // MATRIX3DXD_H