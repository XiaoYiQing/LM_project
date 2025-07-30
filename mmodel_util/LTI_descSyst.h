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




class LTI_descSyst{ 


public:

// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

LTI_descSyst();

LTI_descSyst( Eigen::MatrixXd& E_in, Eigen::MatrixXd& A_in, Eigen::MatrixXd& B_in, 
    Eigen::MatrixXd& C_in, Eigen::MatrixXd& D_in );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Operations
// ====================================================================== >>>>>

/*
Check if the system matrices are consistent with transfer function matrices requirements.
Returns string describing what exactly caused the inconsistency, if any.
*/
string consistency_check() const;
/*
Check if the system matrices are consistent with transfer function matrices requirements.
*/
bool is_consistent() const;

// Obtain the number of outputs.
unsigned int get_output_cnt() const;
// Obtain the number of inputs.
unsigned int get_input_cnt() const;
// Obtain the order of the system.
unsigned int get_order() const;
/* 
Check stability of the system.
WARNING: Depending on the order of the system, this can be costly in computation time.
*/
bool is_stable();
/*
Compute the poles of the system, if possible, and return them.
WARNING: Depending on the order of the system, this can be costly in computation time.
*/
Eigen::VectorXcd get_poles();

/*
Transform the current descriptor system into a regular system (E becomes identity).
WARNING: this operation is irreversible, and is only possible if E is invertible.

Return true if successful.
*/
bool to_reg_syst();

/*
Generate the sparse version of this system, where E becomes the identity matrix and 
A becomes a diagonal matrix (Transfer function is still identical).
WARNING: this operation starts with regular system translation, which is irreversible
if successful.

Return true if successful.
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

WARNING: if the sparse representation has not yet been calculated, this 
fnction call will start the sparse system computation and could add unintended
computation time.
Furthermore, this function fails and return an empty matrix if the system cannot
be sparsified.
*/
Eigen::MatrixXcd tf_sparse_eval( complex<double> );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Function
// ====================================================================== >>>>>

Eigen::MatrixXd get_E() const;
Eigen::MatrixXd get_A() const;
Eigen::MatrixXd get_B() const;
Eigen::MatrixXd get_C() const;
Eigen::MatrixXd get_D() const;

void set_E( const shared_ptr< const Eigen::MatrixXd > E_in );
void set_A( const shared_ptr< const Eigen::MatrixXd > A_in );
void set_B( const shared_ptr< const Eigen::MatrixXd > B_in );
void set_C( const shared_ptr< const Eigen::MatrixXd > C_in );
void set_D( const shared_ptr< const Eigen::MatrixXd > D_in );

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

