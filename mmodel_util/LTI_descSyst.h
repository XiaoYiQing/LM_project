#ifndef LTI_DESCSYST_H
#define LTI_DESCSYST_H


#include <Eigen/Dense>
#include <iostream>
#include <magic_enum.hpp>
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
string consistency_check();

// Obtain the number of outputs.
unsigned int get_output_cnt();
// Obtain the number of inputs.
unsigned int get_input_cnt();
// Obtain the order of the system.
unsigned int get_order();

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Operations
// ====================================================================== >>>>>



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

