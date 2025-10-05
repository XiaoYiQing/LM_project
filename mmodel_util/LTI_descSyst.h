#ifndef LTI_DESCSYST_H
#define LTI_DESCSYST_H

#include <chrono>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <magic_enum.hpp>
#include <memory>
#include <string>   
#include <vector> 

#include "Matrix3DXcd.h"


using namespace std;

// TODO: Put exception throwing instead of simple messages for error handling in this class.

/**
 * Class representing instances of Linear-Time-Invariant (LTI) descriptor systems.
 * It's primary task is holding the E, A, B, C, D matrices and facilitating evaluation
 * of the transfer function on demand.
 */
class LTI_descSyst{


public:

// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

LTI_descSyst();

/**
 * Initialize an LTI_descSyst instance with the given E, A, B, C, D matrices.
 */
LTI_descSyst( Eigen::MatrixXd& E_in, Eigen::MatrixXd& A_in, Eigen::MatrixXd& B_in, 
    Eigen::MatrixXd& C_in, Eigen::MatrixXd& D_in );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Operations
// ====================================================================== >>>>>

/**
 * Check if the system matrices are consistent with transfer function matrices requirements.
 * Returns string describing what exactly caused the inconsistency, if any.
 * 
 * @return String describing the inconsistency. Empty if no inconsistencies.
 */
string consistency_check() const;

/**
 * Obtain the current consistency status.
 * 
 * @return The current consistency flag.
 * 
 * @note Calling this function automatically triggers a consistency check before the flag is returned.
 */
bool is_consistent() const;

/**
 * Obtain the number of outputs (The number of rows to the C matrix).
 * 
 * @return The number of outputs.
 */
unsigned int get_output_cnt() const;
/**
 * Obtain the number of inputs (The number of columns to the B matrix).
 * 
 * @return The number of inputs.
 */
unsigned int get_input_cnt() const;
/**
 * Obtain the order of the system.
 * 
 * @return The order of the system.
 */
unsigned int get_order() const;

/**
 * Check stability of the system.
 * 
 * @return The stability flag of the current system.
 * 
 * @note If the poles aren't currently up to date, this function triggers a pole 
 *  computation process and can be computationally expensive depending on the system order.
 */
bool is_stable();
/**
 * Compute the poles of the system, if possible, and return them.
 * 
 * @return The vector containing the poles of the system.
 * 
 * @note If the poles aren't currently up to date, this function triggers a pole 
 *  computation process and can be computationally expensive depending on the system order.
 */
Eigen::VectorXcd get_poles();

//TODO: make this function not return bool, but throw exceptions if something goes wrong.
/**
 * Transform the current descriptor system into a regular system (E becomes identity).
 * 
 * @return Boolean indicating whether this function has successfully converted the descriptor 
 *  system into a regular system.
 * 
 * @note This operation is irreversible, and is only possible if E is invertible.
 */
bool to_reg_syst();

//TODO: make this function not return bool, but throw exceptions if something goes wrong.
/**
 * Generate the sparse version of this system, where E becomes the identity matrix and 
 * A becomes a diagonal matrix (Transfer function is still identical).
 * 
 * @return Boolean indicating whether the sparsification process is successful.
 * 
 * @note This operation starts with regular system translation, which is irreversible
 *  if successful.
 */
bool gen_sparse_syst();

/**
 * Evaluate the transfer function represented by the current system at the 
 * target frequency.
 * 
 * @param f_tar The frequency evaluation point to applied to the transfer function.
 * @return The evaluated frequency parameter matrix at the evaluation point.
 * 
 * @note No assumption is made on the nature of the frequency value "f_tar", and it will be
 *  used directly in the transfer function equation.
 * 
 * @warning Transfer function evaluation can be expensive depending on the order of the system.
 *  Consider using the sparse transfer function evaluation, if permissible.
 */
Eigen::MatrixXcd tf_eval( complex<double> f_tar ) const;

/**
 * Evaluate the transfer function represented by the current system at the 
 * target frequencies.
 * 
 * @param f_vec Vector of frequency evaluation points to be used in the transfer function.
 * @return Vector of evaluated frequency parameter matrices.
 * 
 * @note No assumption is made on the nature of the frequency values in "f_vec", and 
 *  they will be used directly in the transfer function equation.
 * 
 * @warning Transfer function evaluation can be expensive depending on the order of the system.
 *  Consider using the sparse transfer function evaluation, if permissible.
 */
Matrix3DXcd tf_eval( vector< complex<double> >& f_vec ) const;

/**
 * Evaluate the transfer function using the sparse representation at the 
 * target frequency.
 * 
 * @param f_tar The frequency evaluation point to applied to the transfer function.
 * @return The evaluated frequency parameter matrix at the evaluation point.
 * 
 * @note No assumption is made on the nature of the frequency value "f_tar", and it will be
 *  used directly in the transfer function equation.
 * 
 * @warning If the sparse representation is not up-to-date, this 
 *  function returns an empty result.
 */
Eigen::MatrixXcd tf_sparse_eval( complex<double> f_tar ) const;

/**
 * Evaluate the transfer function the sparse representation at the 
 * target frequencies.
 * 
 * @param f_vec Vector of frequency evaluation points to be used in the transfer function.
 * @return Vector of evaluated frequency parameter matrices.
 * 
 * @note No assumption is made on the nature of the frequency values in "f_vec", and 
 *  they will be used directly in the transfer function equation.
 * 
 * @warning If the sparse representation is not up-to-date, this 
 *  function returns an empty result.
 */
Matrix3DXcd tf_sparse_eval( vector< complex<double> >& f_vec ) const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Function
// ====================================================================== >>>>>

// Return the current E matrix.
Eigen::MatrixXd get_E() const;
// Return the current A matrix.
Eigen::MatrixXd get_A() const;
// Return the current B matrix.
Eigen::MatrixXd get_B() const;
// Return the current C matrix.
Eigen::MatrixXd get_C() const;
// Return the current D matrix.
Eigen::MatrixXd get_D() const;

/**
 * @brief Set the E matrix.
 * @param E_in The new E matrix to be applied.
 * @note This function invalidates up-to-date flags for computed poles and sparse system.
 */
void set_E( const Eigen::MatrixXd& E_in );
/**
 * @brief Set the A matrix.
 * @param E_in The new A matrix to be applied.
 * @note This function invalidates up-to-date flags for computed poles and sparse system.
 */
void set_A( const Eigen::MatrixXd& A_in );
/**
 * @brief Set the B matrix.
 * @param E_in The new B matrix to be applied.
 * @note This function invalidates up-to-date flags for computed poles and sparse system.
 */
void set_B( const Eigen::MatrixXd& B_in );
/**
 * @brief Set the C matrix.
 * @param E_in The new C matrix to be applied.
 * @note This function invalidates up-to-date flags for computed poles and sparse system.
 */
void set_C( const Eigen::MatrixXd& C_in );
/**
 * @brief Set the D matrix.
 * @param E_in The new D matrix to be applied.
 * @note This function invalidates up-to-date flags for computed poles and sparse system.
 */
void set_D( const Eigen::MatrixXd& D_in );

// Obtain the up-to-date flag of computed poles (true if current poles are up-to-date).
bool get_utd_poles() const;
// Obtain the up-to-date flag of sparse system (true if current sparse system is up-to-date).
bool get_utd_sparse_syst() const;

// TODO: You need to prevent the following functions to work if sparse system is not utd.
/**
 * Obtain the diagonalized A matrix.
 * 
 * @return The diagonalized A matrix.
 */
Eigen::SparseMatrix< complex<double> > get_As() const;
/**
 * Obtain the left sparse system transformation matrix.
 * 
 * @return The left sparse system transformation matrix.
 */
Eigen::MatrixXcd get_Ts_L() const;
/**
 * Obtain the right sparse system transformation matrix.
 * 
 * @return The right sparse system transformation matrix.
 */
Eigen::MatrixXcd get_Ts_R() const;

// ====================================================================== <<<<<

protected:

// ====================================================================== >>>>>
//      Member Variables
// ====================================================================== >>>>>

    // Matrix E.
    Eigen::MatrixXd E;
    // Matrix A.
    Eigen::MatrixXd A;
    // Matrix B.
    Eigen::MatrixXd B;
    // Matrix C.
    Eigen::MatrixXd C;
    // Matrix D.
    Eigen::MatrixXd D;

    

    // Boolean indicating if the current computed poles are up-to-date.
    bool utd_poles;
    // Vector holding the poles of the system. Initially empty, but is filled
    // when appropriate function is called.
    Eigen::VectorXcd poles;

    // Boolean indicating if the sparse version is up-to-date.
    bool utd_sparse_syst;
    // The sparse matrix A.
    Eigen::SparseMatrix< complex<double> > As;
    // The left transformation matrix for sparsifying the transfer function.
    Eigen::MatrixXcd Ts_L;
    // The right transformation matrix for sparsifying the transfer function.
    Eigen::MatrixXcd Ts_R;


// ====================================================================== <<<<<

};





#endif  // LTI_DESCSYST_H

