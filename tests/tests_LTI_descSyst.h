#ifndef TESTS_LTI_DESCSYST
#define TESTS_LTI_DESCSYST

#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <string>

#include "LTI_descSyst.h"
#include "numUtils.h"

extern string RES_PATH_XYQ_str;


namespace tests{


    /*
    Test base functionalities.
    case_idx:
        0: base functionalities.
        1: stability.
        2:
    */
    void LTI_descSyst_test_1( unsigned int case_idx );

    /*
    Test regular system translation and sparsification.
    */
    void LTI_descSyst_test_2( unsigned int case_idx );



}



#endif  // TESTS_LTI_DESCSYST