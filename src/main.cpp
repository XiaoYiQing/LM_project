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
#include "tests_Eigen.h"
#include "tests_fData.h"
#include "tests_LM_eng.h"
#include "tests_LTI_descSyst.h"
#include "tests_Matrix3DXd.h"
#include "tests_numUtils.h"



using namespace std;

// Global string variable holding the full path of the extra resource directory.
string RES_PATH_XYQ_str = string( RES_PATH_XYQ );
// Global string variable holding the full path of the source directory.
string SRC_PATH_XYQ_str = string( SRC_PATH_XYQ );


void singVal_extract_run();


int main() {

    int err_ret_val = 1;

    // Test aspectes of the Eigen library.
    // int testCase = 3;

    // tests::eigen_test1( testCase );
    // tests::eigen_test2( 1 );


    // tests::numUtils_test_1(4);
    // tests::eigenUtils_test_1(0);

    // tests::Matrix3DXd_test_1(1);
    // tests::Matrix3DXd_test_2_ops(3);
    // tests::Matrix3DXd_test_3_spec_ops(0);
    // tests::Matrix3DXd_test_4_supp(0);

    // tests::fData_test_1( 0 );
    // tests::fData_test_2( 3 );
    // tests::fData_enum_test(0);
    // tests::fData_test_sXp_read(1);
    // tests::fData_setFunc_tests( 1 );
    // tests::fData_LTspice_data_read_test();
    // tests::fData_print_test(0);

    // tests::LM_eng_test_1( 0 );
    // tests::LM_eng_test_2( 0 );
    // tests::LM_eng_test_3( 1 );
    // tests::LM_eng_class_test();
    
    // tests::LTI_descSyst_test_1( 3 );
    // tests::LTI_descSyst_test_2( 2 );

    // tests::LM_eng_full_SFML_testrun_v2();
    // tests::LM_eng_full_SFML_testrun_gen();
    // tests::LM_eng_full_SFML_dc_case_run();
    tests::LM_eng_print_singVals();

    // singVal_extract_run();

    return 0; 

}





void singVal_extract_run(){

    // RES_PATH_XYQ_str

    vector<string> file_stem_arr;
    vector<string> file_path_arr;
    vector<string> file_ext_arr;
    vector<fData::FDATA_FORMAT> fData_format_arr;
    vector<fData::METRIC_PREFIX> fData_metPrefix_arr;

    unsigned int file_cnt = 0;

    file_stem_arr.push_back( "Slink_a=100um_b=400um" );
    file_path_arr.push_back( RES_PATH_XYQ_str + "/" );
    file_ext_arr.push_back( ".s2p" );
    fData_format_arr.push_back( fData::FDATA_FORMAT::RI );
    fData_metPrefix_arr.push_back( fData::METRIC_PREFIX::G );
    file_cnt++;

    file_stem_arr.push_back( "Slink_a=100um_b=425um" );
    file_path_arr.push_back( RES_PATH_XYQ_str + "/" );
    file_ext_arr.push_back( ".s2p" );
    fData_format_arr.push_back( fData::FDATA_FORMAT::RI );
    fData_metPrefix_arr.push_back( fData::METRIC_PREFIX::G );
    file_cnt++;

    string fullFileName_z = "";
    for( unsigned int z = 0; z < file_cnt; z++ ){

        // Assemble current full filename.
        fullFileName_z = file_path_arr.at(z) + file_stem_arr.at(z) + file_ext_arr.at(z);

        // Define our frequency data object.
        fData myFData;
        // Obtain the data from the target data file and insert into the fData object.
        fData::read_sXp_file( myFData, fullFileName_z );

    }


    int lol = 0;

}
