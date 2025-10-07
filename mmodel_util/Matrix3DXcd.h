#ifndef MATRIX3DXCD_H
#define MATRIX3DXCD_H


#include <Eigen/Dense>
#include <iostream>
#include <magic_enum.hpp>
#include <stdexcept>
#include <string>   
#include <vector> 

#include "Matrix3DXd.h"

using namespace std;



/**
 * A pseudo extension the the Eigen::Matrix2DXcd that is not official.
 * This class mimics a 3D matrix, but is in actuality a vector of 2D matrices 
 * (vector of Eigen::Matrix2DXcd).
 * More than just a vector of 2D matrices, this class allows basic matrix operations
 * seemlessly, which is the primary goal of this class in the first place.
 * The 3rd dimension being represented as the vector dimension also allows fairly easy
 * modification on the 3rd dimension.
 */
class Matrix3DXcd{    


public:

// ====================================================================== >>>>>
//      Static Functions
// ====================================================================== >>>>>

    static const double DEF_NUM_THRESH;


    /**
     * Check the consistency of a Matrix3DXcd instance.
     * For example, all 2D matrices in the vector must have the same size.
     * 
     * @param tarMat The target Matrix3DXcd instance to check for consistency.
     * @return Consistency boolean of the target Matrix3DXcd instance.
     */
    static bool consist_check( const Matrix3DXcd& tarMat );
    /**
     * Check the consistency of a vector of 2D matrices (Eigen::MatrixXcd).
     * For example, all 2D matrices in the vector must have the same size.
     * 
     * @param tarMat The target vector of MatrixXcd to check for consistency.
     * @return Consistency boolean of the target vector of MatrixXcd.
     */
    static bool consist_check( const vector< Eigen::MatrixXcd >& );

    /**
     * Support function for checking whether two Eigen::MatrixXcd instances
     * have the same dimensions.
     * 
     * @param matA 2D matrix number 1.
     * @param matB 2D matrix number 2.
     * @return Boolean indicating whether two two matrices have the same dimensions.
     */
    static bool same_size( const Eigen::MatrixXcd& matA, const Eigen::MatrixXcd& matB );


    /**
     * Function checks whether the 3D matrix contains 0 row and 0 column matrices, in
     * which case the Matrix3DXcd is considered null.
     * 
     * @param tarMat The target matrix to check for null status.
     * @return Boolean indicating whether the Matrix3DXcd instance should be considered null.
     */
    static bool null_ref_check( const Matrix3DXcd& tarMat );
    /*
    Function for checking if the reference matrix of the matrix vector has
    0 rows or 0 columns.
    If true, the reference matrix has 0 rows or columns.
    */
    /**
     * Function checks whether the vector of MatrixXcd contains 0 row and 0 column MatrixXcd.
     * 
     * @param tarMat The target MatrixXcd vector to check for null status.
     * @return Boolean indicating whether the MatrixXcd vector should be considered null.
     */
    static bool null_ref_check( const vector< Eigen::MatrixXcd >& );


// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

    Matrix3DXcd();

    /**
     * Construct a Matrix3DXcd instance as a vector zeros matrices with specified 
     * dimensions.
     * 
     * @param row_idx Number of rows.
     * @param col_idx Number of columns.
     * @param lvl_idx Number of 2D matrices or height or depth of the 3D matrix.
     */
    Matrix3DXcd( unsigned int row_idx, unsigned int col_idx, unsigned int lvl_idx );

    /**
     * Construct a Matrix3DXcd by directly assigning the vector of MatrixXcd.
     * 
     * @param Mat3D The vector of MatrixXcd to be directly integrated into the Matrix3DXcd instance.
     */
    Matrix3DXcd( const vector< Eigen::MatrixXcd >& Mat3D );

    /**
     * Construct a Matrix3DXcd by directly using one Matrix3DXd instance as the real part and
     * one Matrix3DXd instance as the imaginary part.
     */
    Matrix3DXcd( const Matrix3DXd& rePart, const Matrix3DXd& imPart );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Operators
// ====================================================================== >>>>>

    /**
     * Addition operator overload: current Matrix3DXcd instance addition with another 
     * Matrix3DXcd instance.
     * 
     * @param tarMat Target Matrix3DXcd being added to current instance.
     * @return Matrix3DXcd sum.
     */
    Matrix3DXcd operator+(const Matrix3DXcd& tarMat) const;

    // Overloading the substraction operator with another matrix of the same class.
    /**
     * Substraction operator overload: another Matrix3DXcd instance substracted from current 
     * Matrix3DXcd instance.
     * 
     * @param tarMat Target Matrix3DXcd being substracted from the current instance.
     * @return Matrix3DXcd difference.
     */
    Matrix3DXcd operator-(const Matrix3DXcd& tarMat) const;

    /**
     * Scalar multiplication operator overload: A scalar multiplied to the current 
     * Matrix3DXcd instance.
     * 
     * @param scalar Target scalar being multiplied to the current instance.
     * @return Matrix3DXcd product.
     */
    Matrix3DXcd operator*(const double scalar) const;

    /**
     * Scalar compound multiplication operator overload: A scalar multiplied to the 
     * current Matrix3DXcd instance and directly updated in said instance.
     * 
     * @param scalar Target scalar being multiplied to the current instance.
     * @return The same Matrix3DXcd instance, but updated as the product.
     */
    Matrix3DXcd& operator*=(const double scalar);

    /**
     * @brief Same class multiplication operator overload: another Matrix3DXcd instance 
     * multiplied with the current instance element-wise.
     * 
     * In actuality, element-wise multiplication is conducted, where each scalar 
     * element in either matrices sharing the same coordinates are multiplied 
     * together and the product is placed at the same coordinate in the product
     * matrix.
     * 
     * @param tarMat Target Matrix3DXcd instance being multiplied to the current instance.
     * @return Matrix3DXcd element-wise product.
     */
    Matrix3DXcd operator*(const Matrix3DXcd& tarMat) const;

    /**
     * @brief Same class division operator overload: current Matrix3DXcd instance divided
     * by another instance element-wise. 
     * 
     * In actuality, element-wise division is conducted, where each scalar element 
     * in the current matrix instance is divided by the corresponding scalar element 
     * in the divisor matrix of the same coordinates. The quotient scalar is placed 
     * in the quotient matrix at the same coordinate.
     * 
     * @param tarMat Target Matrix3DXcd acting as the divisor to the current isntance.
     * @return Matrix3DXcd element-wise quotient.
     */
    Matrix3DXcd operator/(const Matrix3DXcd& tarMat) const;

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
    void elem_log10();

    /*
    Perform element-wise cosine function on the entries of the matrix.
    */
    void elem_cos();

    /*
    Perform element-wise cosine function on the entries of the matrix.
    */
    void elem_sin();

    /*
    Perform element-wise arctan function on the entries of the matrix.
    */
    void elem_atan();

    /*
    Perform element-wise division with entries of "tarMat" serving as divisors.
    This function has a special rule where if both the entry of the present matrix
    and the one in tarMat at the matching coordinate are zero (according to num_thresh), 
    the result is 0 and not considered a runtime error.
    Context: this is to deal with DC point frequency data having zero real part.
    */
    Matrix3DXcd elem_div_spec( const Matrix3DXcd& tarMat ) const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Operations
// ====================================================================== >>>>>

    /*
    Specialize function which computes the phase (in radians) with the individual
    vector real parts in "rePart" and corresponding imaginary parts in "imPart".
    Each entry in the 3D matrix is considered one independent vector.
    The resulting phases are stored in a 3D matrix of the same size.
    */
    // static Matrix3DXcd elem_phase_comp( const Matrix3DXcd& rePart, const Matrix3DXcd& imPart );
    static Matrix3DXcd elem_phase_comp( const Matrix3DXcd& tarMat );

    /*
    Compute the root-mean-square (RMS) value from the 3D matrix of complex value
    with real and imaginary parts stored separately.
    */
    // static double RMS_total_comp( const Matrix3DXcd& rePart, const Matrix3DXcd& imPart );
    static double RMS_total_comp( const Matrix3DXcd& tarMat );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Support Functions
// ====================================================================== >>>>>

    /* 
    Check if the target 2D matrix can be added to the 3D matrix.
    */
    bool mat3DValidInCheck( const Eigen::MatrixXcd& val );

    /*
    Re-initialize the 3D matrix with an initiate set of zero matrices of specified sizes.
    This function effectively serves as the "reserve" function specifically when the 3D matrix
    is empty (or emptied) and requires a new start with specificed number of rows and columns.
    WARNING: This function erases all current 3D matrix entries.
    */
    void reInit( unsigned int row_cnt, unsigned int col_cnt, unsigned int lvl_cnt );

    /*
    Clear all existing 2D matrices.
    */
    void clear();

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

    /*
    Obtain a continuous subset segment of the vector of 2D matrices.
    */
    Matrix3DXcd segment( size_t startIndex, size_t len ) const;

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
    bool isEmpty() const;

    /*
    Return the 2D matrix at the target index of the vector of 2D matrix.
    NOTE: this "at( unsigned int )" function does not perform the assignment 
    function because it doesn't return a reference, but a copy. This is to 
    prevent insertion of a 2D matrix having different dimensions than the rest.
    */
    Eigen::MatrixXcd at( unsigned int ) const;

    /*
    Return a subset 3D matrix comprised of the group of 2D matrices given by the
    input index vector.
    */
    Matrix3DXcd at( vector< unsigned int > idxVec );

    /*
    Set the 2D matrix values at the target index to the given 2D matrix values.
    */
    void set( unsigned int, const Eigen::MatrixXcd& );

    /*
    Set the value at the specific coordinate of the specific 2D matrix. 
    */
    void set( unsigned int row_idx, unsigned int col_idx, unsigned int lvl_idx, double val );

    /*
    Put a new 2D matrix entry at the end of the vector.
    */
    void push_back( const Eigen::MatrixXcd& in2DMat );
    /*
    Put a vector of 2D matrix entries at the end of the vector.
    */
    void push_back( const vector< Eigen::MatrixXcd >& in2DMatVec );

    /*
    Insert a new 2D matrix entry at the target index position.
    */
    void insert( unsigned int idx, const Eigen::MatrixXcd& in2DMat );


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
    vector< Eigen::MatrixXcd > Mat3D;

// ====================================================================== <<<<<

};



#endif  // MATRIX3DXCD_H