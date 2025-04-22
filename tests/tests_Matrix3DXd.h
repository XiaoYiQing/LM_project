#ifndef TESTS_MATRIX3DXD_H
#define TESTS_MATRIX3DXD_H


#include <iostream>
#include <string>   
#include <vector> 

#include "Matrix3DXd.h"



using namespace std;


namespace tests{

    /*
    Tests focused on basic functionalities.
    */
    void Matrix3DXd_test_1( int case_idx );

    /*
    Tests focused on operators and operations.
    */
    void Matrix3DXd_test_2_ops( int case_idx );

    /*
    Tests focused on specialize operators.
    */
    void Matrix3DXd_test_2_spec_ops( int case_idx );


}







#endif  // TESTS_MATRIX3DXD_H
