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

}

#endif  // TESTS_LM_ENG_H