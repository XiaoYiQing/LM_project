#ifndef LTI_DESCSYST_H
#define LTI_DESCSYST_H


#include <Eigen/Dense>
#include <iostream>
#include <magic_enum.hpp>
#include <memory>
#include <string>   
#include <vector> 



using namespace std;




class LTI_descSyst{ 


public:

// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

LTI_descSyst();

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Operations
// ====================================================================== >>>>>

/*
Check if the system matrices are consistent with transfer function matrices requirements.
*/
string consistency_check() const;

// Obtain the number of outputs.
unsigned int get_output_cnt() const;
// Obtain the number of inputs.
unsigned int get_input_cnt() const;
// Obtain the order of the system.
unsigned int get_order() const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Function
// ====================================================================== >>>>>

Eigen::MatrixXd get_E() const;
Eigen::MatrixXd get_A() const;
Eigen::MatrixXd get_B() const;
Eigen::MatrixXd get_C() const;
Eigen::MatrixXd get_D() const;

Eigen::MatrixXd set_E( const shared_ptr< const Eigen::MatrixXd > E_in );
Eigen::MatrixXd set_A();
Eigen::MatrixXd set_B();
Eigen::MatrixXd set_C();
Eigen::MatrixXd set_D();

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

// ====================================================================== <<<<<

};





#endif  // LTI_DESCSYST_H

