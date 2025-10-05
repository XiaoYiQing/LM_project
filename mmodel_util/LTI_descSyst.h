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

/**
 * Transform the current descriptor system into a regular system (E becomes identity).
 * 
 * @return Boolean indicating whether this function has successfully converted the descriptor 
 *  system into a regular system.
 * 
 * @note This operation is irreversible, and is only possible if E is invertible.
 */
bool to_reg_syst();

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

/*
Evaluate the transfer function represented by the current system at the 
target frequency.

WARNING: transfer function evaluation can be expensive depending on the order of the system.
    Consider using the sparse transfer function evaluation, if permissible.
*/
Eigen::MatrixXcd tf_eval( complex<double> ) const;

/*
Evaluate the transfer function represented by the current system at the 
target frequency.
Input is vector of complex<double>.
Output is array of MatrixXcd matrices, under the class Matrix3DXcd.

WARNING: transfer function evaluation can be expensive depending on the order of the system.
    Consider using the sparse transfer function evaluation, if permissible.
*/
Matrix3DXcd tf_eval( vector< complex<double> >& ) const;

/*
Evaluate the transfer function using the sparse representation at the target
frequency.

WARNING: if the sparse representation is not up-to-date, this 
function returns an empty result.
*/
Eigen::MatrixXcd tf_sparse_eval( complex<double> ) const;

/*
Evaluate the transfer function represented by the current system at the 
target frequency.
Input is vector of complex<double>.
Output is array of MatrixXcd matrices, under the class Matrix3DXcd.

WARNING: if the sparse representation is not up-to-date, this 
function returns an empty result.
*/
Matrix3DXcd tf_sparse_eval( vector< complex<double> >& ) const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Function
// ====================================================================== >>>>>

Eigen::MatrixXd get_E() const;
Eigen::MatrixXd get_A() const;
Eigen::MatrixXd get_B() const;
Eigen::MatrixXd get_C() const;
Eigen::MatrixXd get_D() const;

void set_E( const Eigen::MatrixXd& E_in );
void set_A( const Eigen::MatrixXd& A_in );
void set_B( const Eigen::MatrixXd& B_in );
void set_C( const Eigen::MatrixXd& C_in );
void set_D( const Eigen::MatrixXd& D_in );

bool get_utd_poles() const;
Eigen::VectorXcd get_poles() const;

bool get_utd_sparse_syst() const;

Eigen::SparseMatrix< complex<double> > get_As() const;
Eigen::MatrixXcd get_Ts_L() const;
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

