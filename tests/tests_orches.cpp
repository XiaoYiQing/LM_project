#include "tests_orches.h"






void tests::SFLM_full_run_test( unsigned int test_idx ){

    // Initialize test case index.
    int case_cnt = 0;

    int case_idx = 0;
    // 0- Simple standard run.
    if( case_cnt == case_idx ){

// ---------------------------------------------------------------------- >>>>>
//      Initialization (Data)
// ---------------------------------------------------------------------- >>>>>
    
        // Define our frequency data object.
        fData myFData;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/test_res_dir/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myFData, fullFileName );

        // Switch the data format into real + imaginary format.
        myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
        // Normalize the frequency vector (As much as you can according to metric prefixes).
        myFData.data_prefix_switch( fData::METRIC_PREFIX::M );

// ---------------------------------------------------------------------- <<<<<

        // Create a custom reduced f set index vector.
        unsigned int test_fr_size = LM_eng::STD_RED_FSET_SIZE + 12;
        vector<unsigned int> f_r_idx_vec = 
            utils::gen_lin_idx_arr( 0, myFData.get_f_cnt() - 1, min( test_fr_size, myFData.get_f_cnt() ) );

        // Generate the interleaving relative partition index arrays.
        vector< vector< unsigned int > > index_arrs = 
            LM_UTIL::gen_2_partit_idx_arr( test_fr_size );
        vector< unsigned int > f1_idx_vec = index_arrs.at(0);
        vector< unsigned int > f2_idx_vec = index_arrs.at(1);

        // Perform the full LM engine process.
        shared_ptr<LM_eng> my_LM_eng = 
            FCT_SCR::SFLM_full_run( myFData, f_r_idx_vec, f1_idx_vec, f2_idx_vec );


        // Generate a LIT system.
        shared_ptr<LTI_descSyst> myTF = my_LM_eng->step5_LM_to_tf( 49 );


// ---------------------------------------------------------------------- >>>>>
//      Evaluation
// ---------------------------------------------------------------------- >>>>>

        // Generate the sparse system equivalent.
        myTF->gen_sparse_syst();

        // Generate a frequency evaluation array.
        Eigen::VectorXd tmp_fvec = myFData.getF_vec();
        vector< complex<double> > testFVec = vector< complex<double> >( tmp_fvec.size() );
        for( int z = 0; z < tmp_fvec.size(); z++ ){
            testFVec[z] = complex<double>( 0.0, tmp_fvec(z) );
        }

        // Evaluate the transfer function over the specified frequency array values.
        Matrix3DXcd H_app_mat_arr = myTF->tf_sparse_eval( testFVec );
        // Obtain the original data as a array of complex matrices.
        Matrix3DXcd H_orig_mat_arr = Matrix3DXcd( myFData.getXr_vec(), myFData.getXi_vec() );
        // Compute the difference between the original and approximated frequency data.
        Matrix3DXcd H_diff = H_orig_mat_arr - H_app_mat_arr;

        // Compute the RMS error.
        double total_RMS_err = Matrix3DXcd::RMS_total_comp( H_diff );
        // cout << "The total RMS error: " << total_RMS_err << endl;
        // Test for stability.
        bool is_stab = myTF->is_stable();
        // cout << "System stability: " << is_stab << endl;

        bool test_bool = true;
        test_bool = test_bool && ( abs( total_RMS_err - 0.00050943459036 ) < 1e-8 );
        test_bool = test_bool && is_stab;

        if( test_bool ){
            cout << "SFLM_full_run test: passed!" << endl;
        }else{
            cout << "SFLM_full_run test: failed!" << endl;
        }


// ---------------------------------------------------------------------- <<<<<

    }

}
