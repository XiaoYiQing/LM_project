#include "tests_LM_eng.h"



using namespace std;

void tests::LM_eng_test_1( unsigned int test_idx ){

    int case_cnt = 0;

// ---------------------------------------------------------------------- >>>>>
//      All Tests Common Initialization 
// ---------------------------------------------------------------------- >>>>>
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
// ---------------------------------------------------------------------- <<<<<

    // 0- LM construct test.
    if( test_idx == case_cnt ){

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

    }


    case_cnt++;
    // 1- sLM construct test.
    if( test_idx == case_cnt ){

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


void tests::LM_eng_test_2( unsigned int test_idx ){

    int case_cnt = 0;

// ---------------------------------------------------------------------- >>>>>
//      All Tests Common Initialization 
// ---------------------------------------------------------------------- >>>>>

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
    // Construct the Loewner Matrix using the two partitions.
    shared_ptr<Eigen::MatrixXcd> mySLM = LM_UTIL::build_SLM( *partit1, *partit2 );

// ---------------------------------------------------------------------- <<<<<

    // 0- Base LM pencil verification.
    if( test_idx == case_cnt ){

        complex<double> cplx_f_ref = partit1->get_cplx_f_at( partit1->get_f_cnt()/2 );
        shared_ptr<Eigen::MatrixXcd> myLM_pen = LM_UTIL::build_LM_pencil( cplx_f_ref, *myLM, *mySLM );

        // Perform SVD.
        Eigen::JacobiSVD<Eigen::MatrixXcd> mySVD( *myLM_pen );
        // Get the singular values
        Eigen::VectorXd singularValues = mySVD.singularValues();
        
        // std::cout << std::fixed << std::setprecision(12);
        // cout << singularValues << endl;

        Eigen::VectorXd top10singVals(10);
        top10singVals << 1482.643133463571, 1468.219478553086, 1418.868597959851,
        1385.038432483787, 1311.652883867499, 1294.188802726333, 1230.153395906074,
        1216.890355807787, 1187.052516244472, 1170.663508169054;

        bool match_bool = true;
        for( unsigned int z = 0; z < 10; z++ ){
            match_bool = match_bool && ( top10singVals(z) - singularValues(z) < 1e-9 );
        }
        cout << "Top singular values match: " << match_bool << endl;

    }


    // 1- 
    if( test_idx == case_cnt ){


    }

}



void tests::LM_eng_test_3( unsigned int test_idx ){

    int case_cnt = 0;

    // 0- Real transformation matrix vector mult check.
    if( test_idx == case_cnt ){

        bool has_DC_pt = false;
        unsigned int sub_mat_size = 3;
        unsigned int sub_blk_cnt = 2;

        shared_ptr<Eigen::MatrixXcd> myTMat = LM_UTIL::build_reT_mat( has_DC_pt, sub_mat_size, sub_blk_cnt );

        // Initialize our test vector matrix.
        Eigen::MatrixXcd myMatVet_A( sub_mat_size, 2*sub_blk_cnt*sub_mat_size );
        Eigen::MatrixXcd myMatVet_B( 2*sub_blk_cnt*sub_mat_size, sub_mat_size );

        // Fill the matrix vector with random entries, but following the real transform structure.
        for( unsigned int z = 0; z < sub_blk_cnt; z++ ){

            // Generate random values for the real and imaginary parts.
            shared_ptr<vector<double>> reVec_A = utils::rDoubleGen( -1, 1, sub_mat_size*sub_mat_size );
            shared_ptr<vector<double>> imVec_A = utils::rDoubleGen( -1, 1, sub_mat_size*sub_mat_size );
            shared_ptr<vector<double>> reVec_B = utils::rDoubleGen( -1, 1, sub_mat_size*sub_mat_size );
            shared_ptr<vector<double>> imVec_B = utils::rDoubleGen( -1, 1, sub_mat_size*sub_mat_size );

            // Create the current sub-matrix.
            Eigen::MatrixXcd mat_A_z( sub_mat_size, sub_mat_size );
            Eigen::MatrixXcd mat_B_z( sub_mat_size, sub_mat_size );
            for( unsigned int y = 0; y < reVec_A->size(); y++ ){
                mat_A_z(y) = complex<double>( reVec_A->at(y), imVec_A->at(y) );
                mat_B_z(y) = complex<double>( reVec_B->at(y), imVec_B->at(y) );
            }

            unsigned int lead_orig = z*2*sub_mat_size;
            unsigned int lead_conj = z*2*sub_mat_size + sub_mat_size;

            myMatVet_A.block( 0, lead_orig, sub_mat_size, sub_mat_size ) = mat_A_z;
            myMatVet_A.block( 0, lead_conj, sub_mat_size, sub_mat_size ) = mat_A_z.conjugate();
            myMatVet_B.block( lead_orig, 0, sub_mat_size, sub_mat_size ) = mat_B_z;
            myMatVet_B.block( lead_conj, 0, sub_mat_size, sub_mat_size ) = mat_B_z.conjugate();

        }

        // Multiply the target matrix with the real transform matrix from the RHS.
        Eigen::MatrixXcd myReMatVet_A = myMatVet_A*( *myTMat );
        Eigen::MatrixXcd myReMatVet_B = ( myTMat->conjugate().transpose() )*myMatVet_B;

        bool match_bool = true;
        match_bool = match_bool && ( myReMatVet_A.imag().cwiseAbs().maxCoeff() < 1e-12 );
        cout << "Real transform matrix check (RHS): " << match_bool << endl;
        match_bool = match_bool && ( myReMatVet_B.imag().cwiseAbs().maxCoeff() < 1e-12 );
        cout << "Real transform matrix check (LHS): " << match_bool << endl;
        

    }

    case_cnt++;
    // 1- Real transformation matrix LM mult check.
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

        unsigned int sub_blk_cnt = 8;
        // Create a subset linear index array.
        vector< unsigned int > fr_idx_arr = utils::gen_lin_idx_arr( 0, myFData.get_f_cnt() - 1, sub_blk_cnt );
        // Create a fData subset.
        shared_ptr<fData> myFr = myFData.red_partit( fr_idx_arr );

        // Generate two partitions from this data subset.
        vector< shared_ptr<fData> > myPartits = myFr->gen_2_partit();
        // Generate the two partitions with their complex conjugates inserted 
        // in interleaving fashion.
        shared_ptr<fData> myFr1 = myPartits.at(0)->gen_cplx_conj_comb();
        shared_ptr<fData> myFr2 = myPartits.at(1)->gen_cplx_conj_comb();

        // Construct the Loewner Matrix using the two cconj injected partitions.
        Eigen::MatrixXcd myLM = *LM_UTIL::build_LM( *myFr1, *myFr2 );
        // Construct the Loewner Matrix using the two cconj injected partitions.
        Eigen::MatrixXcd mySLM = *LM_UTIL::build_SLM( *myFr1, *myFr2 );

        // Obtain base parameters of the two partitions.
        bool f1_has_DC_pt = myFr1->hasDC();
        bool f2_has_DC_pt = myFr1->hasDC();
        unsigned int sub_mat_size = myFr1->get_out_cnt();
        // Partition sizes before cconj injection.
        unsigned int sub_blk_cnt_1 = myPartits.at(0)->get_f_cnt();  
        unsigned int sub_blk_cnt_2 = myPartits.at(1)->get_f_cnt();

        // Build the left and right transformation matrices.
        Eigen::MatrixXcd myTMat_L = *LM_UTIL::build_reT_mat( f2_has_DC_pt, sub_mat_size, sub_blk_cnt_1 );
        Eigen::MatrixXcd myTMat_R = *LM_UTIL::build_reT_mat( f1_has_DC_pt, sub_mat_size, sub_blk_cnt_2 );
        // Obtain the hermitian of the right transform matrix.
        Eigen::MatrixXcd myTMat_R_herm = myTMat_R.conjugate().transpose();

        // Perform the transformation.
        Eigen::MatrixXcd myLM_re = ( myTMat_R_herm*myLM )*myTMat_L;


        bool match_bool = true;
        match_bool = match_bool && ( myLM_re.imag().cwiseAbs().maxCoeff() < 1e-12 );
        cout << "Real transform matrix full matrix mult check: " << match_bool << endl;

    }

}



void tests::LM_eng_full_SFML_testrun(){

    

}
