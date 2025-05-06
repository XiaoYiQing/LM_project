#ifndef TESTS_LM_ENG_H
#define TESTS_LM_ENG_H


#include <Eigen/Dense>
#include <iostream>
#include <string>

#include "fData.h"
#include "LM_eng.h"


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
    Special test which performs the entire SFML (System Format Loewner Matrix) process
    and check for mistakes.
    */
    void LM_eng_full_SFML_run();

}

#endif  // TESTS_LM_ENG_H