// Preprocessor variable created when running CMAKE process. This variable 
// indicates the location where non-compilabe or non c/c++ support files are located.
#ifndef RES_PATH_XYQ
#   define RES_PATH_XYQ ""
#else
#endif

// Preprocessor variable created when running CMAKE process which indicates
// the source directory.
#ifndef SRC_PATH_XYQ
#   define SRC_PATH_XYQ ""
#else
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <filesystem>
#include <fstream>
#include "Matrix3DxD.h"
#include <regex>
#include <sstream>
#include <string>   
#include <vector>  



#include "fData.h"
#include "numUtils.h"
#include "orchestrators.h"
#include "tests_eigenUtils.h"
#include "tests_Eigen.h"
#include "tests_fData.h"
#include "tests_LM_eng.h"
#include "tests_LTI_descSyst.h"
#include "tests_Matrix3DXd.h"
#include "tests_numUtils.h"
#include "tests_orches.h"



using namespace std;

// Global string variable holding the full path of the extra resource directory.
string RES_PATH_XYQ_str = string( RES_PATH_XYQ );
// Global string variable holding the full path of the source directory.
string SRC_PATH_XYQ_str = string( SRC_PATH_XYQ );



int main() {

    int err_ret_val = 1;

    // Test aspectes of the Eigen library.
    // int testCase = 3;

    // tests::eigen_test1( 3 );
    // tests::eigen_test2( 1 );
    
    // tests::gen_rand_MatrixXd_test();
    // tests::gen_randn_MatrixXd_test();
    // tests::gen_orth_basis_test();
    // tests::SVD_econ_test();
    // tests::rSVD_test();

    // tests::Vd_to_file_test();
    // tests::file_to_vec_test();
    // tests::MatrixXd_to_file_test();
    // tests::file_to_MatrixXd_test();
    // tests::MatrixXcd_to_file_test();
    // tests::file_to_MatrixXcd_test();

    // tests::numUtils_test_1(4);
    // tests::gen_match_vector_test();
    // tests::sort_num_vec_inplace_test();

    // tests::Matrix3DXd_test_1(5);
    // tests::Matrix3DXd_test_2_ops(3);
    // tests::Matrix3DXd_test_3_spec_ops(0);
    // tests::Matrix3DXd_test_4_supp(0);
    // tests::Matrix3DXd_test_serialize();

    // tests::fData_test_1( 0 );
    // tests::fData_test_2( 3 );
    // tests::fData_enum_test(0);
    // tests::fData_test_sXp_read(1);
    // tests::fData_setFunc_tests( 1 );
    // tests::fData_LTspice_data_read_test();
    // tests::fData_print_test(0);
    // tests::fData_prefix_manip_test(1);
    // tests::fData_serialize_test();

    // tests::LM_eng_cplx_LM_test();
    // tests::LM_pencil_test();
    // tests::LM_eng_reT_test();
    // tests::LM_eng_serialize_test();
    
    tests::LTI_descSyst_test_1( 1 );
    // tests::LTI_descSyst_access_consist_test();
    // tests::LTI_descSyst_test_2( 2 );

    // tests::LM_eng_full_SFML_testrun();
    // tests::LM_eng_full_SFML_testrun_v2();
    // tests::LM_eng_full_SFML_testrun_gen();
    // tests::LM_eng_class_test();
    // tests::LM_eng_full_SFML_dc_case_run();
    // tests::LM_eng_rSVD_case_run();
    // tests::LM_eng_rSVD_case_run_vb();
    // tests::LM_eng_print_singVals();
    // tests::LM_eng_re_LM_comp_test( 1 );

    // tests::orch_SFLM_full_run_test();
    // tests::orch_SFLM_direc_re_run_test();
    // FCT_SCR::singVal_extract_run();

    return 0; 

}




