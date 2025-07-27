#include "LTI_descSyst.h"




// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>


LTI_descSyst::LTI_descSyst(){

    /*
    Create the default transfer function matrices that are nothing but
    matrices of size 1x1 and containing just 1.
    */
    this->E = Eigen::MatrixXd::Ones(1,1);
    this->A = Eigen::MatrixXd::Ones(1,1);
    this->B = Eigen::MatrixXd::Ones(1,1);
    this->C = Eigen::MatrixXd::Ones(1,1);
    this->D = Eigen::MatrixXd::Zero(1,1);

}

// ====================================================================== <<<<<




// ====================================================================== >>>>>
//      Specialized Operations
// ====================================================================== >>>>>

bool LTI_descSyst::is_consistent() const{

    string tmp_str = this->consistency_check();
    return tmp_str.empty();

}

/*
Check if the system matrices are consistent with transfer function matrices requirements.
A consistent system return an empty string.
*/
string LTI_descSyst::consistency_check() const{

    string err_msg = "";
    if( this->A.rows() != this->E.rows() || this->A.cols() != this->E.cols() ){
        err_msg = "E and A matrices are not consistent in size.";
        return err_msg;
    }

    if( this->A.rows() != this->B.rows() ){
        err_msg = "Matrix B's size is inconsistent with E and A matrices.";
        return err_msg;
    }

    if( this->A.cols() != this->C.cols() ){
        err_msg = "Matrix C's size is inconsistent with E and A matrices.";
        return err_msg;
    }

    if( this->C.rows() != this->D.rows() || this->B.cols() != this->D.cols() ){
        err_msg = "Matrix D's size is inconsistent with C and B matrices.";
        return err_msg;
    }

    return err_msg;

}

// Obtain the number of outputs.
unsigned int LTI_descSyst::get_output_cnt() const{
    return static_cast<unsigned int>( this->B.cols() );
}
// Obtain the number of inputs.
unsigned int LTI_descSyst::get_input_cnt() const{
    return static_cast<unsigned int>( this->C.rows() );
}
// Obtain the order of the system.
unsigned int LTI_descSyst::get_order() const{
    return static_cast<unsigned int>( this->A.rows() );
}

bool LTI_descSyst::is_stable() const{

    // Verify if the current system is legitimate.
    if( !this->is_consistent() ){
        cout << "System is inconsistent: cannot generate poles." << endl;
        return false;
    }

    // Compute E^(-1)*A as solution x to E*x = A
    auto solver = this->E.fullPivHouseholderQr(); // or other suitable decomposition
    Eigen::MatrixXd pole_mat = solver.solve( this->A );

    // Compute eigenvalues of E^(-1)*A
    Eigen::ComplexEigenSolver< Eigen::MatrixXcd > mySolver( pole_mat );
    // Check if the computation was successful
    if ( mySolver.info() != Eigen::Success ) {
        std::cerr << "Failed to compute eigenvalues." << std::endl;
        return false;
    }

    // Determine if the system is stable (Maximum poles real part is negative).
    Eigen::VectorXcd eigeVals_1 = mySolver.eigenvalues();
    bool is_stab = 0 > eigeVals_1.real().maxCoeff();

    return is_stab;

}


Eigen::MatrixXcd LTI_descSyst::tf_eval( complex<double> f_tar ) const{

    Eigen::MatrixXcd tmp_1 = ( f_tar * this->E - this->A );
    Eigen::MatrixXcd B_tmp = this->B;

    // Solve system using LU decomposition
    Eigen::FullPivLU<Eigen::MatrixXcd> lu( tmp_1 );
    if(!lu.isInvertible()) {
        std::cerr << "Matrix is singular or nearly singular." << std::endl;
        return -1*Eigen::MatrixXcd::Ones(1,1);
    }

    Eigen::MatrixXcd tmp_2 = lu.solve( B_tmp );
    Eigen::MatrixXcd H_z = this->C * tmp_2;
    
    return H_z;

}


Matrix3DXcd LTI_descSyst::tf_eval( vector< complex<double> > f_vec ) const{

    unsigned int eval_cnt = f_vec.size();

    Matrix3DXcd res_var = Matrix3DXcd( this->get_input_cnt(), this->get_output_cnt(),
        eval_cnt );

    Eigen::MatrixXcd B_tmp = this->B;

    for( unsigned int z = 0; z < eval_cnt; z++ ){

        Eigen::MatrixXcd tmp_1 = ( f_vec[z] * this->E - this->A );

        // Solve system using LU decomposition
        Eigen::FullPivLU<Eigen::MatrixXcd> lu( tmp_1 );
        if(!lu.isInvertible()) {
            std::cerr << "Matrix is singular or nearly singular." << std::endl;
            return Matrix3DXcd(1,1,1);
        }
        Eigen::MatrixXcd tmp_2 = lu.solve( B_tmp );
        Eigen::MatrixXcd H_z = this->C * tmp_2;
        
        res_var.set( z, H_z );

    }

    return res_var;

}

// ====================================================================== <<<<<



// ====================================================================== >>>>>
//      Access Function
// ====================================================================== >>>>>

Eigen::MatrixXd LTI_descSyst::get_E() const{
    return this->E;
}
Eigen::MatrixXd LTI_descSyst::get_A() const{
    return this->A;
}
Eigen::MatrixXd LTI_descSyst::get_B() const{
    return this->B;
}
Eigen::MatrixXd LTI_descSyst::get_C() const{
    return this->C;
}
Eigen::MatrixXd LTI_descSyst::get_D() const{
    return this->D;
}

void LTI_descSyst::set_E( const shared_ptr< const Eigen::MatrixXd > E_in )
    { this->E = *E_in; }
void LTI_descSyst::set_A( const shared_ptr< const Eigen::MatrixXd > A_in )
    { this->A = *A_in; }
void LTI_descSyst::set_B( const shared_ptr< const Eigen::MatrixXd > B_in )
    { this->B = *B_in; }
void LTI_descSyst::set_C( const shared_ptr< const Eigen::MatrixXd > C_in )
    { this->C = *C_in; }
void LTI_descSyst::set_D( const shared_ptr< const Eigen::MatrixXd > D_in )
    { this->D = *D_in; }

// ====================================================================== <<<<<
