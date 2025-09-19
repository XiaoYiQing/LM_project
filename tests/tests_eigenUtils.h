#ifndef TEST_EIGENUTILS_H
#define TEST_EIGENUTILS_H


#include <chrono>
#include <iostream>
#include <string>

#include "eigenUtils.h"
#include "numUtils.h"


using namespace std;


extern string RES_PATH_XYQ_str;
extern string SRC_PATH_XYQ_str;

namespace tests{


    // Test functions in the eigenUtils file.
    void eigenUtils_test_1( unsigned int test_idx );

    void Vd_to_file_test();

    /*
    Test the function file_to_vec()
    */
    void file_to_vec_test();

    void MatrixXd_to_file_test();

    void file_to_MatrixXd_test();

    void MatrixXcd_to_file_test();

    void file_to_MatrixXcd_test();

    void gen_rand_MatrixXd_test();
    
    void gen_randn_MatrixXd_test();

    void gen_orth_basis_test();

    void SVD_econ_test();

}


#endif  // TEST_EIGENUTILS_H