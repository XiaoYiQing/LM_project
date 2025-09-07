#include "tests_eigenUtils.h"


void tests::eigenUtils_test_1( unsigned int test_idx ){

    // Initialize test case index.
    int case_cnt = 0;

    int case_idx = 0;
    // 0- Test gen_lin_idx_arr.
    if( case_cnt == case_idx ){

        unsigned int num_cnt = 10;

        shared_ptr< vector<double> > tmp_vec = utils::rDoubleGen( -10.0, 10.0, num_cnt );

        // Map the std::vector data to Eigen::VectorXd
        Eigen::Map<Eigen::VectorXd> eigenVec(tmp_vec->data(), tmp_vec->size());

        // If you need an independent copy
        Eigen::VectorXd tarVec = eigenVec;

        string targetDir = "C:/Users/Yi Qing Xiao/Documents/Cpp_projects/LM_project/data_output";
        string targetStemName = "tmp_data_file";

        utils::vec_to_file( targetDir, targetStemName, tarVec, 0 );

    }

}

