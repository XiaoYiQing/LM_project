#include "tests_LM_eng.h"



using namespace std;

void tests::LM_eng_test_1( unsigned int test_idx ){

    int case_cnt = 0;

    // 0- LM construct test.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myFData;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myFData, fullFileName );

        // Switch the data format into real + imaginary format.
        myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
        // Normalize the frequency vector (As much as you can according to metric prefixes).
        myFData.data_prefix_switch( fData::METRIC_PREFIX::G );

        // Create a subset linear index array.
        vector< unsigned int > fr_idx_arr = utils::gen_lin_idx_arr( 0, myFData.get_f_cnt() - 1, 100 );
        // Create a fData subset.
        shared_ptr<fData> myFr = myFData.red_partit( fr_idx_arr );

        // Generate two partitions from this data subset.
        vector< shared_ptr<fData> > myPartits = myFr->gen_2_partit();
        shared_ptr<fData> partit1 = myPartits.at(0);
        shared_ptr<fData> partit2 = myPartits.at(1);
        // Construct the Loewner Matrix using the two partitions.
        shared_ptr<Eigen::MatrixXcd> myLM = LM_UTIL::build_LM( *partit1, *partit2 );

        // Obtain the number of outputs and inputs.
        unsigned int out_cnt = myFr->get_out_cnt();
        unsigned int in_cnt = myFr->get_in_cnt();

        bool match_bool = true;
        // Random test point selection.
        unsigned int test_j = randIntGen( 0, partit1->get_f_cnt() - 1, 1 )->at(0);
        unsigned int test_i = randIntGen( 0, partit2->get_f_cnt() - 1, 1 )->at(0); 
        Eigen::MatrixXcd test_S_j = partit1->get_cplxData_at_f( test_j );
        Eigen::MatrixXcd test_S_i = partit2->get_cplxData_at_f( test_i );
        complex<double> test_f_j = partit1->get_cplx_f_at( test_j );
        complex<double> test_f_i = partit2->get_cplx_f_at( test_i );

        // Compute the current reference LM sub-block.
        Eigen::MatrixXcd test_LM_ij = ( test_S_i - test_S_j )/( test_f_i - test_f_j );
        // Obtain the generated LM sub-block.
        Eigen::MatrixXcd LM_ij = 
            myLM->block( ( test_i )*out_cnt, ( test_j )*in_cnt, out_cnt, in_cnt );
        // Compute difference between the two LM sub-blocks.
        Eigen::MatrixXcd LM_ij_diff = test_LM_ij - LM_ij;

        match_bool = match_bool && ( LM_ij_diff.cwiseAbs().maxCoeff() < 1e-9 );

        cout << "Random LM sub-block (" << test_i << ", " << test_j << ") match: " 
            << match_bool << endl;

        // // Compute the SVD
        // Eigen::JacobiSVD<Eigen::MatrixXcd> mySVD( myLM, Eigen::ComputeThinU | Eigen::ComputeThinV );
        // // Get the singular values
        // Eigen::VectorXd singularValues = mySVD.singularValues();
        // cout << singularValues << endl;

    }


    case_cnt++;
    // 1- sLM construct test.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myFData;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myFData, fullFileName );

        // Switch the data format into real + imaginary format.
        myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
        // Normalize the frequency vector (As much as you can according to metric prefixes).
        myFData.data_prefix_switch( fData::METRIC_PREFIX::G );

        // Create a subset linear index array.
        vector< unsigned int > fr_idx_arr = utils::gen_lin_idx_arr( 0, myFData.get_f_cnt() - 1, 100 );
        // Create a fData subset.
        shared_ptr<fData> myFr = myFData.red_partit( fr_idx_arr );

        // Generate two partitions from this data subset.
        vector< shared_ptr<fData> > myPartits = myFr->gen_2_partit();
        shared_ptr<fData> partit1 = myPartits.at(0);
        shared_ptr<fData> partit2 = myPartits.at(1);
        // Construct the shifted Loewner Matrix using the two partitions.
        shared_ptr<Eigen::MatrixXcd> mySLM = LM_UTIL::build_SLM( *partit1, *partit2 );

        // Obtain the number of outputs and inputs.
        unsigned int out_cnt = myFr->get_out_cnt();
        unsigned int in_cnt = myFr->get_in_cnt();

        bool match_bool = true;
        // Random test point selection.
        unsigned int test_j = randIntGen( 0, partit1->get_f_cnt() - 1, 1 )->at(0);
        unsigned int test_i = randIntGen( 0, partit2->get_f_cnt() - 1, 1 )->at(0); 
        Eigen::MatrixXcd test_S_j = partit1->get_cplxData_at_f( test_j );
        Eigen::MatrixXcd test_S_i = partit2->get_cplxData_at_f( test_i );
        complex<double> test_f_j = partit1->get_cplx_f_at( test_j );
        complex<double> test_f_i = partit2->get_cplx_f_at( test_i );

        // Compute the current reference LM sub-block.
        Eigen::MatrixXcd test_SLM_ij = ( test_f_i*test_S_i - test_f_j*test_S_j )/( test_f_i - test_f_j );
        // Obtain the generated LM sub-block.
        Eigen::MatrixXcd SLM_ij = 
            mySLM->block( ( test_i )*out_cnt, ( test_j )*in_cnt, out_cnt, in_cnt );
        // Compute difference between the two LM sub-blocks.
        Eigen::MatrixXcd SLM_ij_diff = test_SLM_ij - SLM_ij;

        match_bool = match_bool && ( SLM_ij_diff.cwiseAbs().maxCoeff() < 1e-9 );

        cout << "Random SLM sub-block (" << test_i << ", " << test_j << ") match: " 
            << match_bool << endl;

    }


    case_cnt++;
    // 2- W construct test.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myFData;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myFData, fullFileName );

        // Switch the data format into real + imaginary format.
        myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
        // Normalize the frequency vector (As much as you can according to metric prefixes).
        myFData.data_prefix_switch( fData::METRIC_PREFIX::G );

        // Create a subset linear index array.
        vector< unsigned int > fr_idx_arr = utils::gen_lin_idx_arr( 0, myFData.get_f_cnt() - 1, 100 );
        // Create a fData subset.
        shared_ptr<fData> myFr = myFData.red_partit( fr_idx_arr );

        // Generate two partitions from this data subset.
        vector< shared_ptr<fData> > myPartits = myFr->gen_2_partit();
        shared_ptr<fData> partit1 = myPartits.at(0);
        shared_ptr<fData> partit2 = myPartits.at(1);
        // Construct the W matrix vector using partition 1.
        shared_ptr<Eigen::MatrixXcd> myW = LM_UTIL::build_W( *partit1 );

        // Obtain the number of outputs and inputs.
        unsigned int out_cnt = myFr->get_out_cnt();
        unsigned int in_cnt = myFr->get_in_cnt();

        bool match_bool = true;
        // Random test point selection.
        unsigned int test_j = randIntGen( 0, partit1->get_f_cnt() - 1, 1 )->at(0);
        Eigen::MatrixXcd test_S_j = partit1->get_cplxData_at_f( test_j );

        // Obtain the generated LM sub-block.
        Eigen::MatrixXcd W_ij = 
            myW->block( 0, ( test_j )*in_cnt, out_cnt, in_cnt );
        // Compute difference between the two LM sub-blocks.
        Eigen::MatrixXcd W_ij_diff = test_S_j - W_ij;

        match_bool = match_bool && ( W_ij_diff.cwiseAbs().maxCoeff() < 1e-9 );

        cout << "Random W sub-block (" << test_j << ") match: " 
            << match_bool << endl;

    }


    case_cnt++;
    // 3- F construct test.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myFData;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myFData, fullFileName );

        // Switch the data format into real + imaginary format.
        myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
        // Normalize the frequency vector (As much as you can according to metric prefixes).
        myFData.data_prefix_switch( fData::METRIC_PREFIX::G );

        // Create a subset linear index array.
        vector< unsigned int > fr_idx_arr = utils::gen_lin_idx_arr( 0, myFData.get_f_cnt() - 1, 100 );
        // Create a fData subset.
        shared_ptr<fData> myFr = myFData.red_partit( fr_idx_arr );

        // Generate two partitions from this data subset.
        vector< shared_ptr<fData> > myPartits = myFr->gen_2_partit();
        shared_ptr<fData> partit1 = myPartits.at(0);
        shared_ptr<fData> partit2 = myPartits.at(1);
        // Construct the F matrix vector using partition 2.
        shared_ptr<Eigen::MatrixXcd> myF = LM_UTIL::build_F( *partit2 );

        // Obtain the number of outputs and inputs.
        unsigned int out_cnt = myFr->get_out_cnt();
        unsigned int in_cnt = myFr->get_in_cnt();

        bool match_bool = true;
        // Random test point selection.
        unsigned int test_i = randIntGen( 0, partit2->get_f_cnt() - 1, 1 )->at(0);
        Eigen::MatrixXcd test_S_i = partit2->get_cplxData_at_f( test_i );

        // Obtain the generated LM sub-block.
        Eigen::MatrixXcd F_ij = 
            myF->block( ( test_i )*out_cnt, 0, out_cnt, in_cnt );
        // Compute difference between the two LM sub-blocks.
        Eigen::MatrixXcd F_ij_diff = test_S_i - F_ij;

        match_bool = match_bool && ( F_ij_diff.cwiseAbs().maxCoeff() < 1e-9 );

        cout << "Random F sub-block (" << test_i << ") match: " 
            << match_bool << endl;

    }

}

void LM_eng_test_2( unsigned int test_idx ){
    
}
