#ifndef ORCHESTRATORS_H
#define ORCHESTRATORS_H



#include <string>
#include <vector>

#include "fData.h"
#include "LM_eng.h"


using namespace std;

extern string RES_PATH_XYQ_str;
extern string SRC_PATH_XYQ_str;


namespace FCT_SCR{

    /*
    Function follows a predefined list of S-param data files and creates a
    LM_eng for each file and then print the real LM pencil singular values to
    a predefined data output directory.
    */
    void singVal_extract_run();

    /*
    Function performs a full SFLM run based on the given data.
    Inputs:
    - src_data: the frequency data object with which the system is going to be 
        constructed and tested.
    - f_r_idx_vec: the index vector with which a subset fData of "src_data" is created.
        This subset fData is the one actually used to construct the LM system.
        The remaining data in "src_data" are used for validation.
    - f1_idx_vec: the index vector creating partition #1 from the subset fData.
    - f2_idx_vec: the index vector creating partition #2 from the subset fData.
    - systOrd: the final transfer function order. Note that this specified number should
        be less than the size of the Loewner Matrix, as your system cannot be bigger
        than what the Loewner Matrix can provide.
    */
    LM_eng SFLM_full_run( const fData& src_data, 
        vector<unsigned int> f_r_idx_vec, vector<unsigned int> f1_idx_vec,
        vector<unsigned int> f2_idx_vec );

    /*
    Function performs a SFLM run skipping computation of complex LMs based on the given 
    data. This means the final LM_eng holds no complex Loewner Matrices.
    Inputs:
    - src_data: the frequency data object with which the system is going to be 
        constructed and tested.
    - f_r_idx_vec: the index vector with which a subset fData of "src_data" is created.
        This subset fData is the one actually used to construct the LM system.
        The remaining data in "src_data" are used for validation.
    - f1_idx_vec: the index vector creating partition #1 from the subset fData.
    - f2_idx_vec: the index vector creating partition #2 from the subset fData.
    - systOrd: the final transfer function order. Note that this specified number should
        be less than the size of the Loewner Matrix, as your system cannot be bigger
        than what the Loewner Matrix can provide.
    */
    LM_eng SFLM_direc_re_run( const fData& src_data, 
        vector<unsigned int> f_r_idx_vec, vector<unsigned int> f1_idx_vec,
        vector<unsigned int> f2_idx_vec );

};



#endif  // ORCHESTRATORS_H