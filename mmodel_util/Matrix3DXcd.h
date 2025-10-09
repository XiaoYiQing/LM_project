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
 * \brief A pseudo unofficial extension the the Eigen::Matrix2DXcd.
 * 
 * This class mimics a 3D matrix, but is in actuality a vector of 2D matrices 
 * (vector of Eigen::Matrix2DXcd).
 * 
 * More than just a vector of 2D matrices, this class allows basic matrix operations
 * seemlessly, which is the primary goal of this class in the first place.
 * 
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
     * @param row_cnt Number of rows.
     * @param col_cnt Number of columns.
     * @param lvl_cnt Number of 2D matrices or height or depth of the 3D matrix.
     */
    Matrix3DXcd( unsigned int row_cnt, unsigned int col_cnt, unsigned int lvl_cnt );

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
     * @param tarMat Target Matrix3DXcd acting as the divisor to the current instance.
     * @return Matrix3DXcd element-wise quotient.
     */
    Matrix3DXcd operator/(const Matrix3DXcd& tarMat) const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Operations
// ====================================================================== >>>>>

    /**
     * Perform element-wise power to the target scalar on all entries of the matrix.
     * 
     * @param exp_val Target exponent to raise all entries of the matrix to the power of.
     */
    void elem_pow( double exp_val );

    /**
     * Perform element-wise power with the input base value raised to the power value given
     * by the entries of the matrix.
     * 
     * @param base_val Target base value to be raied to the power of each entries of the matrix.
     */
    void elem_raise_pow( double base_val );

    /**
     * Perform element-wise base 10 logarithmic operation.
     */
    void elem_log10();

    /**
     * Perform element-wise cosine function on the entries of the matrix. 
     * 
     * @note Entries considered in radians.
     */
    void elem_cos();

    /**
     * Perform element-wise sine function on the entries of the matrix. 
     * 
     * @note Entries considered in radians.
     */
    void elem_sin();

    /**
     * Perform element-wise arctan function on the entries of the matrix. 
     * 
     * @note Entries considered in radians.
     */
    void elem_atan();
    
    /**
     * @brief Perform element-wise division with entries of "tarMat" serving as 
     * divisors with special rules regarding division by 0.
     * 
     * This function has a special rule where if both the entry of the present 
     * matrix and the one in tarMat at the matching coordinate are zero 
     * (according to num_thresh), the result is 0 and not considered a runtime
     * error.
     * 
     * @param tarMat Target Matrix3DXcd acting as the divisor to the current instance.
     * @return Matrix3DXcd element-wise quotient.
     * 
     * @note This function is designed to deal with DC point frequency data having zero real part.
     */
    Matrix3DXcd elem_div_spec( const Matrix3DXcd& tarMat ) const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Operations
// ====================================================================== >>>>>

    /**
     * @brief Specialized function which computes the phase (in radians) of each 
     * element of the target matrix.
     * 
     * @param tarMat Target matrix for which the element-wise phases are to be computed.
     * @return The 3D matrix containing the element-wise phases.
     */
    static Matrix3DXcd elem_phase_comp( const Matrix3DXcd& tarMat );

    /**
     * Compute the root-mean-square (RMS) value from the target Matrix3DXcd instance.
     * 
     * @param tarMat Target matrix from which the RMS is to be computed.
     * @return The RMS value of the target matrix.
     */
    static double RMS_total_comp( const Matrix3DXcd& tarMat );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Support Functions
// ====================================================================== >>>>>

    /**
     * Check if the target 2D matrix can be added to the 3D matrix.
     * 
     * @param val Target 2D matrix for which entry into the current 3D matrix instance is checked.
     * @return Boolean indicating whether the target 2D matrix can be added or inserted in the
     *  current instance.
     */
    bool mat3DValidInCheck( const Eigen::MatrixXcd& val );

    /**
     * @brief Re-initialize the 3D matrix with a set of zero matrices of specified sizes.
     * 
     * This function effectively serves as the "reserve" function specifically 
     * when the 3D matrix is empty (or emptied) and requires a new start with 
     * specificed number of rows and columns, though this function differs in 
     * that actual memory are dedicated when calling it.
     * 
     * @param row_cnt Number of rows.
     * @param col_cnt Number of columns.
     * @param lvl_cnt Number of 2D matrices or height or depth of the 3D matrix.
     * 
     * @warning This function erases all current 3D matrix entries.
     * 
     */
    void reInit( unsigned int row_cnt, unsigned int col_cnt, unsigned int lvl_cnt );

    /**
     * @brief Clear all existing 2D matrices.
     * 
     * This is achieved by emptying the vector holding the 2D matrices.
     */
    void clear();


    /**
     * \brief Reserve a set number of entries worth of memory for additional 2D matrices.
     * 
     * Unlike the real "reserve" function which allocates memory without creating matrix entries,
     * this function creates a number of zeros matrices in the vector to fulfill the requested
     * total number of 2D matrices to reserve.
     * 
     * \param res_size The total amount of 2D matrices the vector of 2D matrices is specified 
     *  to hold.
     * 
     * \note res_size is not the number of additional entries to add to the current 2D matrix vector;
     *  res_size is the final size of the 2D matrix vector. If the current vector has 6 entries and 
     *  res_size = 10, then the vector is going to be given 4 more (zero) matrix entries.
     *  If res_size <= current 2D matrix vector size, nothing happens.
     */
    void reserve( unsigned int res_size );

    /**
     * \brief Set the number of entries of the matrix to the specified size.
     * 
     * Effectively, we simply use the "resize()" function of a vector in general.
     * This function is more definitive than the "reserve" function as it will 
     * delete 2D matrix entries beyond the new vector limit to fulfill the new 
     * size requested.
     * 
     * @param tarSize The amount of entries the 2D matrix vector is forced to contain.
     */
    void resize( unsigned int tarSize );

    /**
     * Perform the same operations as shrink_to_fit() of a vector variable.
     * 
     * @note With how my version of "reserve" functions, this function will basically 
     *  do nothing because I am not actually using the true "reserve" function of a vector.
     */
    void shrink_to_fit();

    /**
     * Obtain a continuous subset segment of the vector of 2D matrices.
     * 
     * @param startIndex The starting index of the subset of the 2D matrix vector to extract.
     * @param len The number of entries to extract for subset.
     * @return The subset version of the original Matrix3DXcd.
     * 
     */
    Matrix3DXcd segment( size_t startIndex, size_t len ) const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>

    /**
     * Obtain the size of the current 3D matrix in a vector of size 3.
     * 
     * @return Vector of unsigned int containing the number of rows, columns, and 
     *  depth/level/height.
     */
    vector<unsigned int> size() const;

    // Fetch the number of rows directly.
    unsigned int rows() const;
    // Fetch the number of columns directly.
    unsigned int cols() const;
    // Fetch the number of depth/level/height directly.
    unsigned int levels() const;

    /**
     * \brief Determine if the current Matrix3DXcd instance is empty.
     * 
     * \return Boolean indicating whether the 3D matrix is empty.
     * 
     * \note Only checks if there are entries in the vector. Will return false even if 
     *  the vector contains nothing but 0 x 0 matrices.
     */
    bool isEmpty() const;

    /**
     * \brief Return the 2D matrix (copy) at the target index of the vector of 2D matrix.
     * 
     * \param tarIdx The index of the target 2D matrix to extract from the vector.
     * \return The 2D matrix (copy) at the specified index.
     * 
     * \note This "at( unsigned int tarIdx )" function does not perform the assignment 
     * function because it doesn't return a reference, but a copy. This is to 
     * prevent insertion of a 2D matrix having different dimensions than the rest.
     */
    Eigen::MatrixXcd at( unsigned int tarIdx ) const;

    /**
     * Return a subset Matrix3DXcd comprised of the group of 2D matrices given by the
     * input index vector.
     * 
     * \param idxVec The vector of indices of the 2D matrices to extract to create the 
     * subset Matrix3DXcd.
     * \return The subset Matrix3DXcd containing only the 2D matrices of the specified indices.
     */
    Matrix3DXcd at( vector< unsigned int > idxVec );

    /**
     * Set the 2D matrix values at the target index to the given 2D matrix values.
     * 
     * \param tarIdx Index of the target 2D matrix to replace.
     * \param tarMat The new 2D matrix to be placed at the target index.
     * 
     * \note Input matrix must match expected 2D matrix dimensions.
     */
    void set( unsigned int tarIdx, const Eigen::MatrixXcd& tarMat );

    /**
     * Set the value at the specific coordinate of the specific 2D matrix. 
     * 
     * \param row_idx Target cell row index.
     * \param col_idx Target cell column index.
     * \param lvl_idx Target cell height/depth/level index.
     * \param val Target cell's new value to be assigned.
     */
    void set( unsigned int row_idx, unsigned int col_idx, unsigned int lvl_idx, complex<double> val );

    /**
     * Put a new 2D matrix entry at the end of the vector.
     * 
     * \param in2DMat New matrix entry to be added at the end of the vector of 2D matrices.
     * 
     * \note Dimensions of new matrix entry must match dimensions of current vector entries.
     */
    void push_back( const Eigen::MatrixXcd& in2DMat );

    /**
     * Put a vector of 2D matrix entries at the end of the vector.
     * 
     * \param in2DMatVec Vector of new matrix entries to be added at the end of 
     *  of the present vector 2D matrices.
     * 
     * \note Dimensions of new matrix entries must match dimensions of current vector entries.
     */
    void push_back( const vector< Eigen::MatrixXcd >& in2DMatVec );

    /**
     * Insert a new 2D matrix entry at the target index position.
     * 
     * \param idx Index where the new matrix entry is to be placed.
     * \param in2DMat New matrix entry to be inserted into the vector.
     * 
     * \note Dimensions of new matrix entry must match dimensions of current vector entries.
     */
    void insert( unsigned int idx, const Eigen::MatrixXcd& in2DMat );

// ====================================================================== <<<<<

    // TODO: write the serialize + deserialize functions for this class
    /**
     * Serialize method to save current object state to a binary file.
     * 
     * \param filename The filename of the binary file where the binary data is to be written to.
     *  Make sure the extension is ".bin".
     * 
     * \note Make sure you use a full filename to know exactly where the file ends up.
     */
    void serialize( const std::string& filename ) const;
    /**
     * Serialize method to save current object state to a binary file.
     * 
     * \param ofs An existing binary file output stream through which the binary data is sent.
     */
    void serialize( std::ofstream& ofs ) const;

    /**
     * Deserialize method to imprint current object state with existing state written
     * in a binary file.
     * 
     * \param filename The filename of the binary file where the binary data is 
     *  to be extracted from. Make sure the extension is ".bin".
     */
    void deserialize( const std::string& filename );
    /**
     * Deserialize method to imprint current object state with existing state from 
     * currently opened binary file input stream.
     * 
     * \param ifs An existing binary file input stream through which the binary data is extracted.
     */
    void deserialize( std::ifstream& ifs );

protected:

    
// ====================================================================== >>>>>
//      Member Variables
// ====================================================================== >>>>>

    /**
     * \brief The numerical threshold used by this object to determine if a numerical 
     * value is to be considered trivial.
     * 
     * As such, this value is also used to evaluate if two numerical values are equal.
     */
    double num_thresh;

    /**
     * The vector of 2D matrices that act as a pseudo 3D matrix within this class.
     */
    vector< Eigen::MatrixXcd > Mat3D;

// ====================================================================== <<<<<

};



#endif  // MATRIX3DXCD_H