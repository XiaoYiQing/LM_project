#ifndef TESTS_FDATA_H
#define TESTS_FDATA_H

#include <cstdio>
#include <Eigen/Dense>
#include <iostream>
#include <string>


#include "fData.h"



extern string RES_PATH_XYQ_str;
extern string SRC_PATH_XYQ_str;


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

    /*
    fData class related enums tests.
    */
    void fData_enum_test( unsigned int test_idx );

    /*
    Test the parser for touchstone files.
    */
    void fData_test_sXp_read( unsigned int test_idx );

    /*
    Test the set functions.
    */
    void fData_setFunc_tests( unsigned int test_idx );

    /*
    Test the parser for LTspice data files.
    */
    void fData_LTspice_data_read_test();

    /*
    Test for writing the data to a file.
    */
    void fData_print_test( unsigned int test_idx );

    /*
    Test for fData support functions performing manipulations on prefix or scale.
    */
    void fData_prefix_manip_test( unsigned int test_idx );

    /*
    Test serialization of fData objects.
    */
    void fData_serialize_test();

}


#endif  // TESTS_FDATA_H