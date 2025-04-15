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
*/
class Matrix3DXd{    


public:

// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

    Matrix3DXd();

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
    vector<unsigned int> size();

    unsigned int rows();
    unsigned int cols();
    unsigned int levels();

// ====================================================================== <<<<<


    /*
    The 3D matrix, which is just a vector of 2D matrices.
    */
    vector< Eigen::MatrixXd > Mat3D;

protected:

    

private:




};





#endif  // MATRIX3DXD_H