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

