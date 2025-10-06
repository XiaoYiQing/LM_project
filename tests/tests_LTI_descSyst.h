#ifndef TESTS_LTI_DESCSYST
#define TESTS_LTI_DESCSYST

#include <chrono>
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

    void LTI_descSyst_access_consist_test();

    void LTI_descSyst_stab_check_test();

    void LTI_descSyst_tf_eval_test();

    void LTI_descSyst_poles_test();

    /*
    Test regular system translation and sparsification.
    */
    void LTI_descSyst_test_2( unsigned int case_idx );



}



#endif  // TESTS_LTI_DESCSYST