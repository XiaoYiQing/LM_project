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

    void LTI_descSyst_access_consist_test();

    void LTI_descSyst_stab_check_test();

    void LTI_descSyst_tf_eval_test();

    void LTI_descSyst_poles_test();

    void LTI_descSyst_sparse_test();

    void LTI_descSyst_regSyst_test();

    void LTI_descSyst_sparSyst_test();



}



#endif  // TESTS_LTI_DESCSYST