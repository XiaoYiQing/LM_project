#ifndef TESTS_LM_ENG_H
#define TESTS_LM_ENG_H


#include <Eigen/Dense>
#include <iostream>
#include <string>

#include "fData.h"
#include "LM_eng.h"
#include "LTI_descSyst.h"


extern string RES_PATH_XYQ_str;


namespace tests{

    /*
    Test basic LM construction.
    */
    void LM_eng_test_1( unsigned int test_idx );
    
    /*
    Test LM pencil.
    */
    void LM_eng_test_2( unsigned int test_idx );

    /*
    Test real transform matrices.
    */
    void LM_eng_test_3( unsigned int test_idx );

    /*
    Special test which performs the entire SFLM (System Format Loewner Matrix) process
    and check for mistakes.
    */
    void LM_eng_full_SFML_testrun();

    /*
    Special test which performs the entire SFLM process and check for mistakes.
    v2: Added transfer function class for more modular coding.
    */
    void LM_eng_full_SFML_testrun_v2();

    /*
    General run for testing various examples.
    */
    void LM_eng_full_SFML_testrun_gen();

    /*
    Test the class object of the LM_eng.
    */
    void LM_eng_class_test();

}

#endif  // TESTS_LM_ENG_H