#ifndef TESTS_FDATA_H
#define TESTS_FDATA_H

#include <Eigen/Dense>
#include <iostream>
#include <string>


#include "fData.h"



extern string RES_PATH_XYQ_str;

using namespace std;


namespace tests{

    /*
    Base functionalities check.
    */
    void fData_test_1( unsigned int test_idx );

    /*
    Specialized function checks.
    */
    void fData_test_2( unsigned int test_idx );

}


#endif  // TESTS_FDATA_H