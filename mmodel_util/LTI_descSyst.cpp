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

Eigen::MatrixXd LTI_descSyst::set_E( const shared_ptr< const Eigen::MatrixXd > E_in ){
    this->E = *E_in;
}

// ====================================================================== <<<<<
