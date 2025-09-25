#ifndef TESTS_LM_ENG_H
#define TESTS_LM_ENG_H

#include <chrono>
#include <Eigen/Dense>
#include <iostream>
#include <string>

#include "eigenUtils.h"
#include "fData.h"
#include "LM_eng.h"
#include "LTI_descSyst.h"


extern string RES_PATH_XYQ_str;
extern string SRC_PATH_XYQ_str;


namespace tests{

    /*
    Test basic LM construction.
    */
    void LM_eng_cplx_LM_test();
    
    /*
    Test LM pencil.
    */
    void LM_pencil_test();

    /*
    Test real transform matrices.
    */
    void LM_eng_reT_test();

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

    /*
    Specific run for the SFLM process when the zero frequency data point is present
    in the initial data set.
    */
    void LM_eng_full_SFML_dc_case_run();

    /*
    Special SFLM process which uses random SVD to avoid performing full SVD.
    */
    void LM_eng_rSVD_case_run();
    
    /*
    Testing the SFLM partial SVD process from an existing serialized file.
    */
    void LM_eng_rSVD_case_run_vb();

    /*
    Test the automatic singular value to file printing function.
    */
    void LM_eng_print_singVals();

    /*
    Test the step functions.
    */
    void LM_eng_steps_test( unsigned int test_idx );

    /*
    Test the serialization + deserialization functions.
    */
    void LM_eng_serialize_test();

}

#endif  // TESTS_LM_ENG_H