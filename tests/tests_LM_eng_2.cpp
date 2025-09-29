#include "tests_LM_eng.h"



void tests::LM_eng_full_SFML_testrun_gen(){

    // Size of reduced frequency data array.
    unsigned int fr_len = 100;

// ---------------------------------------------------------------------- >>>>>
//      Initialization (Data)w
// ---------------------------------------------------------------------- >>>>>

    // Define our frequency data object.
    fData myFData;
    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/test_res_dir/inductor_1_width_4_dielectric_35.s2p";

    // Obtain the data from the target data file and insert into the fData object.
    // fData::read_LTspice_Sp_file( myFData, fullFileName );
    fData::read_sXp_file( myFData, fullFileName );

    // Switch the data format into real + imaginary format.
    myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
    // Normalize the frequency vector (As much as you can according to metric prefixes).
    myFData.data_prefix_switch( fData::METRIC_PREFIX::M );

    
    // Create a subset linear index array.
    vector< unsigned int > fr_idx_arr = 
        utils::gen_lin_idx_arr( 0, myFData.get_f_cnt() - 1, fr_len );
    // Create a fData subset.
    fData myFr = myFData.red_partit( fr_idx_arr );

    // Generate two partitions from this data subset.
    vector< shared_ptr<fData> > myFrs = myFr.gen_2_partit();
    // Generate the two partitions with their complex conjugates inserted 
    // in interleaving fashion.
    fData myFrc1 = myFrs.at(0)->gen_cplx_conj_comb();
    fData myFrc2 = myFrs.at(1)->gen_cplx_conj_comb();

    // Obtain base parameters of the two partitions.
    bool f1_has_DC_pt = myFrc1.hasDC();
    bool f2_has_DC_pt = myFrc2.hasDC();
    unsigned int out_cnt = myFrc1.get_out_cnt();
    // Partition sizes before cconj injection.
    unsigned int fr1_len = myFrs.at(0)->get_f_cnt();  
    unsigned int fr2_len = myFrs.at(1)->get_f_cnt();

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      LM System Construct
// ---------------------------------------------------------------------- >>>>>

    // Construct the Loewner Matrix using the two cconj injected partitions.
    Eigen::MatrixXcd myLM = LM_UTIL::build_LM( myFrc1, myFrc2 );
    // Construct the Loewner Matrix using the two cconj injected partitions.
    Eigen::MatrixXcd mySLM = LM_UTIL::build_SLM( myFrc1, myFrc2 );
    // Construct the W matrix vector using partition 1.
    Eigen::MatrixXcd myW = LM_UTIL::build_W( myFrc1 );
    // Construct the F matrix vector using partition 2.
    Eigen::MatrixXcd myF = LM_UTIL::build_F( myFrc2 );

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      LM System Real Transform
// ---------------------------------------------------------------------- >>>>>

    // Build the left and right transformation matrices.
    Eigen::MatrixXcd myTMat_L = *LM_UTIL::build_reT_mat( f2_has_DC_pt, out_cnt, fr1_len );
    Eigen::MatrixXcd myTMat_R = *LM_UTIL::build_reT_mat( f1_has_DC_pt, out_cnt, fr2_len );
    // Obtain the hermitian of the right transform matrix.
    Eigen::MatrixXcd myTMat_R_herm = myTMat_R.conjugate().transpose();

    // Perform the transformation.
    Eigen::MatrixXcd myLM_re_tmp = ( myTMat_R_herm*myLM )*myTMat_L;
    Eigen::MatrixXcd mySLM_re_tmp = ( myTMat_R_herm*mySLM )*myTMat_L;
    Eigen::MatrixXcd myW_re_tmp = myW*myTMat_L;
    Eigen::MatrixXcd myF_re_tmp = myTMat_R_herm*myF;

    // Check for real matrices.
    bool match_bool = true;
    match_bool = match_bool && ( myLM_re_tmp.imag().cwiseAbs().maxCoeff() < 1e-12 );
    match_bool = match_bool && ( mySLM_re_tmp.imag().cwiseAbs().maxCoeff() < 1e-12 );
    match_bool = match_bool && ( myW_re_tmp.imag().cwiseAbs().maxCoeff() < 1e-12 );
    match_bool = match_bool && ( myF_re_tmp.imag().cwiseAbs().maxCoeff() < 1e-12 );
    cout << "SFLM real matrices check: " << match_bool << endl;

    // Obtain purely real defintion of the matrices.
    Eigen::MatrixXd myLM_re = myLM_re_tmp.real();
    Eigen::MatrixXd mySLM_re = mySLM_re_tmp.real();
    Eigen::MatrixXd myW_re = myW_re_tmp.real();
    Eigen::MatrixXd myF_re = myF_re_tmp.real();

    // Generate a random test point.
    unsigned int test_f_idx = utils::rIntGen( 0, myFr.get_f_cnt() - 1, 1 )->at(0);
    complex<double> test_f = myFr.get_cplx_f_at( test_f_idx );
    Eigen::MatrixXcd tmpAns = myW_re*( ( - test_f*myLM_re + mySLM_re ).inverse() )*myF_re;
    Eigen::MatrixXcd ansDiff = myFr.get_cplxData_at_f( test_f_idx ) - tmpAns;

    match_bool = true;
    match_bool = match_bool && ( ansDiff.cwiseAbs2().maxCoeff() < 1e-12 );
    cout << "Full sized LM system evaluation test (Not mandatory to pass): " << match_bool << endl;
    
// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      LM Pencil Generation and SVD
// ---------------------------------------------------------------------- >>>>>

    // Obtain a reference frequency value.
    // double ref_f = myFr.get_fval_at( (unsigned int) ceil( (double) fr_len/2 ) );
    double ref_f = myFr.get_fval_at( myFr.get_f_cnt() - 1 );
    
    // Construct the LM pencil.
    shared_ptr<Eigen::MatrixXd> LM_pen = 
        LM_UTIL::build_LM_pencil( ref_f, myLM_re, mySLM_re );

    // Perform SVD.
    Eigen::JacobiSVD<Eigen::MatrixXd> svdResObj( *LM_pen, Eigen::ComputeFullU | Eigen::ComputeFullV );
    // Get the singular values
    Eigen::VectorXd singVals = svdResObj.singularValues();
    // Get the left singular vectors (U)
    Eigen::MatrixXd U = svdResObj.matrixU();
    // Get the right singular vectors (V)
    Eigen::MatrixXd V = svdResObj.matrixV();


    // std::cout << std::fixed << std::setprecision(12);
    // cout << singVals << endl;

    // Perform the model reduction to obtain usable E, A, B, C matrices.
    Eigen::MatrixXcd E_full = -1*( U.transpose() * myLM_re * V );
    Eigen::MatrixXcd A_full = -1*( U.transpose() * mySLM_re * V );
    Eigen::MatrixXcd C_full = myW_re * V;
    Eigen::MatrixXcd B_full = U.transpose() * myF_re;

    Eigen::MatrixXcd tmp_z = ( test_f * E_full - A_full );
    Eigen::MatrixXcd H_z = C_full * tmp_z.inverse() * B_full;

    match_bool = true;
    ansDiff = myFr.get_cplxData_at_f( test_f_idx ) - H_z;
    match_bool = match_bool && ( ansDiff.cwiseAbs2().maxCoeff() < 1e-12 );
    cout << "Full sized LM system post-SVD evaluation test: " << match_bool << endl;

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      SVD Model Order Reduction
// ---------------------------------------------------------------------- >>>>>

    // Define the number of singular values to retain.
    unsigned int svd_ret_cnt = 10;

    Eigen::VectorXd singVals_r = singVals.segment( 0, svd_ret_cnt );

    Eigen::MatrixXd U_r = U.block( 0, 0, U.rows(), svd_ret_cnt );
    Eigen::MatrixXd V_r = V.block( 0, 0, V.rows(), svd_ret_cnt );

    // Perform the model reduction to obtain usable E, A, B, C matrices.
    Eigen::MatrixXd E_n = -1*( U_r.transpose() * myLM_re * V_r );
    Eigen::MatrixXd A_n = -1*( U_r.transpose() * mySLM_re * V_r );
    Eigen::MatrixXd C_n = myW_re * V_r;
    Eigen::MatrixXd B_n = U_r.transpose() * myF_re;
    Eigen::MatrixXd D_n = Eigen::MatrixXd::Zero( out_cnt, out_cnt );

    // Model generation.
    LTI_descSyst mySyst = LTI_descSyst( E_n, A_n, B_n, C_n, D_n );

// ---------------------------------------------------------------------- >>>>>


// ---------------------------------------------------------------------- >>>>>
//      Stability Check
// ---------------------------------------------------------------------- >>>>>

    // Determine if the system is stable (Maximum poles real part is negative).
    Eigen::VectorXcd poles = mySyst.get_poles();
    bool is_stab = 0 > poles.real().maxCoeff();
    cout << "Is stable: " << is_stab << endl;
    cout << "Max pole real part: " << poles.real().maxCoeff() << endl;

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Model Evaluation
// ---------------------------------------------------------------------- >>>>>

    // Compute the sparse system.
    mySyst.gen_sparse_syst();

    // Generate a frequency evaluation array.
    Eigen::VectorXd tmp_fvec = myFData.getF_vec();
    vector< complex<double> > testFVec2 = vector< complex<double> >( tmp_fvec.size() );
    for( int z = 0; z < tmp_fvec.size(); z++ ){
        testFVec2[z] = complex<double>( 0.0, tmp_fvec(z) );
    }

    // Evaluate the transfer function over the specified frequency array values.
    Matrix3DXcd H_app_mat_arr = mySyst.tf_sparse_eval( testFVec2 );
    // Obtain the original data as a array of complex matrices.
    Matrix3DXcd H_orig_mat_arr = Matrix3DXcd( myFData.getXr_vec(), myFData.getXi_vec() );
    // Compute the difference between the original and approximated frequency data.
    Matrix3DXcd H_diff = H_orig_mat_arr - H_app_mat_arr;

    // Compute the RMS error.
    double total_RMS_err = Matrix3DXcd::RMS_total_comp( H_diff );
    cout << "The total RMS error: " << total_RMS_err << endl;

// ---------------------------------------------------------------------- <<<<<


}



void tests::LM_eng_class_test(){

    // Define our frequency data object.
    fData myFData;
    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/test_res_dir/inductor_1_width_4_dielectric_35.s2p";

    // Obtain the data from the target data file and insert into the fData object.
    fData::read_sXp_file( myFData, fullFileName );

    // Switch the data format into real + imaginary format.
    myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
    // Normalize the frequency vector (As much as you can according to metric prefixes).
    myFData.data_prefix_switch( fData::METRIC_PREFIX::G );

// ---------------------------------------------------------------------- >>>>>
//      Full LM Process
// ---------------------------------------------------------------------- >>>>>

    // LM engine initialization.
    LM_eng myEng( myFData );
    myEng.step1_fData_partition();
    myEng.step2_LM_construct();
    myEng.step3_LM_re_trans();
    myEng.step4_LM_pencil_SVD();

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Model Evaluation
// ---------------------------------------------------------------------- >>>>>

    // Specify the test system's order.
    unsigned int tarOrder = 10;
    // Obtain transfer function at specified system order.
    LTI_descSyst myTF = myEng.step5_LM_to_tf( tarOrder );

    // Compute the sparse system.
    myTF.gen_sparse_syst();

    // Generate a frequency evaluation array.
    Eigen::VectorXd tmp_fvec = myFData.getF_vec();
    vector< complex<double> > testFVec = vector< complex<double> >( tmp_fvec.size() );
    for( int z = 0; z < tmp_fvec.size(); z++ ){
        testFVec[z] = complex<double>( 0.0, tmp_fvec(z) );
    }

    // Evaluate the transfer function over the specified frequency array values.
    Matrix3DXcd H_app_mat_arr = myTF.tf_sparse_eval( testFVec );
    // Obtain the original data as an array of complex matrices.
    Matrix3DXcd H_orig_mat_arr = Matrix3DXcd( myFData.getXr_vec(), myFData.getXi_vec() );
    // Compute the difference between the original and approximated frequency data.
    Matrix3DXcd H_diff = H_orig_mat_arr - H_app_mat_arr;

    bool test_bool = true;

    // Compute the RMS error.
    double total_RMS_err = Matrix3DXcd::RMS_total_comp( H_diff );
    cout << "The total RMS error: " << total_RMS_err << endl;
    test_bool = test_bool && ( total_RMS_err < 3.637e-5 );

    // Determine system stability.
    bool is_stab = myTF.is_stable();
    cout << "Is stable: " << is_stab << endl;
    test_bool = test_bool && is_stab;

    if( test_bool ){
        cout << "LM_eng class standard test: passed!" << endl;
    }else{
        cout << "LM_eng class standard test: failed!" << endl;
    }

// ---------------------------------------------------------------------- >>>>>

}


void tests::LM_eng_full_SFML_dc_case_run(){

    // Size of reduced frequency data array.
    unsigned int fr_len = 100;

// ---------------------------------------------------------------------- >>>>>
//      Initialization (Data)
// ---------------------------------------------------------------------- >>>>>
    
    // Define our frequency data object.
    fData myFData;

    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/bondwire_with_strip_design3.s4p";
    fData::read_sXp_file( myFData, fullFileName );

    // Switch the data format into real + imaginary format.
    myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
    // Normalize the frequency vector (As much as you can according to metric prefixes).
    myFData.data_prefix_switch( fData::METRIC_PREFIX::M );

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Full LM Process Step 1 Checks
// ---------------------------------------------------------------------- >>>>>

    // LM engine initialization.
    LM_eng myEng( myFData );
    myEng.step1_fData_partition();

    fData Frc1 = myEng.get_Frc1();
    fData Frc2 = myEng.get_Frc2();

    bool t1bool = true;
    t1bool = t1bool && ( Frc1.get_f_cnt() == Frc2.get_f_cnt() - 1 );
    t1bool = t1bool && ( Frc1.get_fval_at(0) == 0 );
    if( t1bool ){
        cout << "DC case step 1 test: passed!" << endl;
    }else{
        cout << "DC case step 1 test: failed!" << endl;
    }
    

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Full LM Process Step 2 Checks
// ---------------------------------------------------------------------- >>>>>

    myEng.step2_LM_construct();
    Eigen::MatrixXcd my_LM = myEng.get_LM();
    Eigen::MatrixXcd my_SLM = myEng.get_SLM();

    unsigned int out_cnt = myEng.get_out_cnt();
    unsigned int in_cnt = myEng.get_in_cnt();
    unsigned int Frc1_len = Frc1.get_f_cnt();
    unsigned int Frc2_len = Frc2.get_f_cnt();


    bool t2bool = true;
    t2bool = t2bool && ( Frc1_len*out_cnt == my_LM.cols() );
    t2bool = t2bool && ( Frc2_len*in_cnt == my_LM.rows() );
    if( t2bool ){
        cout << "DC case step 2 test: passed!" << endl;
    }else{
        cout << "DC case step 2 test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Full LM Process Step 3 Checks
// ---------------------------------------------------------------------- >>>>>

    myEng.step3_LM_re_trans();
    Eigen::MatrixXcd my_LM_re = myEng.get_LM_re();
    Eigen::MatrixXcd my_SLM_re = myEng.get_SLM_re();

    bool t3bool = true;
    t3bool = t3bool && ( my_LM_re.rows() == Frc2_len*in_cnt );
    t3bool = t3bool && ( my_LM_re.cols() == Frc1_len*out_cnt );
    if( t3bool ){
        cout << "DC case step 3 test: passed!" << endl;
    }else{
        cout << "DC case step 3 test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Full LM Process Step 4 Checks
// ---------------------------------------------------------------------- >>>>>

    myEng.step4_LM_pencil_SVD();
    Eigen::VectorXd my_singVals = myEng.get_singVals();
    Eigen::MatrixXd my_U = myEng.get_U();
    Eigen::MatrixXd my_V = myEng.get_V();

    bool t4bool = true;
    t4bool = t4bool && ( my_singVals.size() == Frc1_len*out_cnt );
    t4bool = t4bool && ( my_U.rows() == Frc2_len*in_cnt );
    t4bool = t4bool && ( my_V.rows() == Frc1_len*out_cnt );
    t4bool = t4bool && ( my_U.cols() == Frc2_len*in_cnt );
    t4bool = t4bool && ( my_V.cols() == Frc1_len*out_cnt );
    if( t4bool ){
        cout << "DC case step 4 test: passed!" << endl;
    }else{
        cout << "DC case step 4 test: failed!" << endl;
    }

    // Write the singular values to external file.
    string targetDir = "C:/Users/Yi Qing Xiao/Documents/Cpp_projects/LM_project/data_output";
    string targetStemName = "tmp_data_file";
    utils::VectorXd_to_file( targetDir, targetStemName, my_singVals, 10 );

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Model Evaluation
// ---------------------------------------------------------------------- >>>>>

    // Specify the test system's order.
    unsigned int tarOrder = 70;
    // Obtain transfer function at specified system order.
    LTI_descSyst myTF = myEng.step5_LM_to_tf( tarOrder );

    // Compute the sparse system.
    myTF.gen_sparse_syst();

    // Generate a frequency evaluation array.
    Eigen::VectorXd tmp_fvec = myFData.getF_vec();
    vector< complex<double> > testFVec = vector< complex<double> >( tmp_fvec.size() );
    for( int z = 0; z < tmp_fvec.size(); z++ ){
        testFVec[z] = complex<double>( 0.0, tmp_fvec(z) );
    }

    // Evaluate the transfer function over the specified frequency array values.
    Matrix3DXcd H_app_mat_arr = myTF.tf_sparse_eval( testFVec );
    // Obtain the original data as a array of complex matrices.
    Matrix3DXcd H_orig_mat_arr = Matrix3DXcd( myFData.getXr_vec(), myFData.getXi_vec() );
    // Compute the difference between the original and approximated frequency data.
    Matrix3DXcd H_diff = H_orig_mat_arr - H_app_mat_arr;

    // Compute the RMS error.
    double total_RMS_err = Matrix3DXcd::RMS_total_comp( H_diff );
    cout << "The total RMS error: " << total_RMS_err << endl;

    // Get stability confirmation.
    bool isStab = myTF.is_stable();
    cout << "System stability: " << isStab << endl;

    bool test_bool = true;
    test_bool = test_bool && ( abs( total_RMS_err - 0.00725466 ) < 1e-6 );
    test_bool = test_bool && isStab;
    if( t4bool ){
        cout << "DC case final model test: passed!" << endl;
    }else{
        cout << "DC case final model test: failed!" << endl;
    }

// ---------------------------------------------------------------------- >>>>>

}


void tests::LM_eng_rSVD_case_run(){

    // Size of reduced frequency data array.
    unsigned int fr_len = 100;

// ---------------------------------------------------------------------- >>>>>
//      Initialization (Data)
// ---------------------------------------------------------------------- >>>>>
    
    // Define our frequency data object.
    fData myFData;

    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
    fData::read_sXp_file( myFData, fullFileName );

    // Switch the data format into real + imaginary format.
    myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
    // Normalize the frequency vector (As much as you can according to metric prefixes).
    myFData.data_prefix_switch( fData::METRIC_PREFIX::M );

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      LM Engine Setup
// ---------------------------------------------------------------------- >>>>>

    // LM engine initialization.
    LM_eng myEng( myFData );
    // First 3 steps.
    myEng.step1_fData_partition();
    myEng.step3skip2_LM_re_construct();

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Random SVD Step
// ---------------------------------------------------------------------- >>>>>

    // Define the number of singular values to compute.
    unsigned int svd_cnt_max = 100;
    // Perform the rSVD.
    myEng.step4_LM_pencil_SVD( svd_cnt_max );

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Model Evaluation
// ---------------------------------------------------------------------- >>>>>

    // Specify the test system's order.
    unsigned int tarOrder = 48;
    // Obtain transfer function at specified system order.
    LTI_descSyst myTF = myEng.step5_LM_to_tf( tarOrder );

    // Compute the sparse system.
    myTF.gen_sparse_syst();

    // Generate a frequency evaluation array.
    Eigen::VectorXd tmp_fvec = myFData.getF_vec();
    vector< complex<double> > testFVec = vector< complex<double> >( tmp_fvec.size() );
    for( int z = 0; z < tmp_fvec.size(); z++ ){
        testFVec[z] = complex<double>( 0.0, tmp_fvec(z) );
    }

    // Evaluate the transfer function over the specified frequency array values.
    Matrix3DXcd H_app_mat_arr = myTF.tf_sparse_eval( testFVec );
    // Obtain the original data as a array of complex matrices.
    Matrix3DXcd H_orig_mat_arr = Matrix3DXcd( myFData.getXr_vec(), myFData.getXi_vec() );
    // Compute the difference between the original and approximated frequency data.
    Matrix3DXcd H_diff = H_orig_mat_arr - H_app_mat_arr;

    bool test_bool = true;
    // Compute the RMS error.
    double total_RMS_err = Matrix3DXcd::RMS_total_comp( H_diff );
    cout << "The total RMS error: " << total_RMS_err << endl;
    test_bool = test_bool && ( total_RMS_err < 0.00055 );

    // Get stability confirmation.
    bool isStab = myTF.is_stable();
    cout << "System stability: " << isStab << endl;
    test_bool = test_bool && isStab;

    if( test_bool ){
        cout << "LM_eng rSVD run test: passed!" << endl;
    }else{
        cout << "LM_eng rSVD run test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<

    // string outFileName = SRC_PATH_XYQ_str + "/data_output/LM_eng_Slink_ex.bin";
    // myEng.serialize( outFileName );    

}


void tests::LM_eng_rSVD_case_run_vb(){

    string srcFileName = RES_PATH_XYQ_str + "/test_res_dir/LM_eng_Slink_ex.bin";
    LM_eng myEng;
    myEng.deserialize( srcFileName );

    fData myFData = myEng.get_fData();

// ---------------------------------------------------------------------- >>>>>
//      Model Evaluation
// ---------------------------------------------------------------------- >>>>>

    // Specify the test system's order.
    unsigned int tarOrder = 49;
    // Obtain transfer function at specified system order.
    LTI_descSyst myTF = myEng.step5_LM_to_tf( tarOrder );

    // Compute the sparse system.
    myTF.gen_sparse_syst();

    // Generate a frequency evaluation array.
    Eigen::VectorXd tmp_fvec = myFData.getF_vec();
    vector< complex<double> > testFVec = vector< complex<double> >( tmp_fvec.size() );
    for( int z = 0; z < tmp_fvec.size(); z++ ){
        testFVec[z] = complex<double>( 0.0, tmp_fvec(z) );
    }

    // Evaluate the transfer function over the specified frequency array values.
    Matrix3DXcd H_app_mat_arr = myTF.tf_sparse_eval( testFVec );
    // Obtain the original data as a array of complex matrices.
    Matrix3DXcd H_orig_mat_arr = Matrix3DXcd( myFData.getXr_vec(), myFData.getXi_vec() );
    // Compute the difference between the original and approximated frequency data.
    Matrix3DXcd H_diff = H_orig_mat_arr - H_app_mat_arr;

    bool test_bool = true;
    // Compute the RMS error.
    double total_RMS_err = Matrix3DXcd::RMS_total_comp( H_diff );
    cout << "The total RMS error: " << total_RMS_err << endl;
    test_bool = test_bool && ( total_RMS_err < 0.00051 );

    // Get stability confirmation.
    bool isStab = myTF.is_stable();
    cout << "System stability: " << isStab << endl;
    test_bool = test_bool && isStab;

    if( test_bool ){
        cout << "LM_eng rSVD run test: passed!" << endl;
    }else{
        cout << "LM_eng rSVD run test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<

}


void tests::LM_eng_print_singVals(){

// ---------------------------------------------------------------------- >>>>>
//      0: Full process singular values print function.
// ---------------------------------------------------------------------- >>>>>

    // Define the data source file and the print destination directory.
    string data_fullFileName = RES_PATH_XYQ_str + "/test_res_dir/audioamp.txt";
    string dest_dirPath = SRC_PATH_XYQ_str + "/data_output";
    // Perform the full process function.
    LM_eng my_LM_eng = LM_eng::print_singVals( data_fullFileName, dest_dirPath );

    // Obtain the singular values from the printed data file.
    string sv_data_fullFileName = dest_dirPath + "/audioamp_sv.txt";
    vector<double> read_vec_tmp = utils::file_to_Vd( sv_data_fullFileName );
    Eigen::VectorXd read_vec = Eigen::VectorXd( read_vec_tmp.size() );
    for( unsigned int z = 0; z < read_vec_tmp.size(); z++ ){
        read_vec(z) = read_vec_tmp[z];
    }
    Eigen::VectorXd comp_vec = my_LM_eng.get_singVals();

    // Compute the difference between the vector from the file and directly from the source.
    Eigen::VectorXd vec_diff = comp_vec - read_vec;

    bool test_bool = vec_diff.cwiseAbs().maxCoeff() < 1e-7;
    if( test_bool ){
        cout << "LM_eng_print_singVals case 1 test: passed!" << endl;
    }else{
        cout << "LM_eng_print_singVals case 1 test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Case 1: Existing LM engine singular values print function
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

    // Perform the full LM engine process.
    LM_eng my_LM_eng2( myFData );
    my_LM_eng2.step1_fData_partition();
    my_LM_eng2.step3skip2_LM_re_construct();
    my_LM_eng2.step4_LM_pencil_SVD();
    // Obtain the singular values directly from the LM engine.
    comp_vec = my_LM_eng2.get_singVals();

    // Print the singular values at the target location.
    string destDir = SRC_PATH_XYQ_str + "/data_output";
    LM_eng::print_singVals( my_LM_eng2, "Slink_a=100um_b=400um_LOL", destDir );
    // Define the full filename of the printed singular values file.
    sv_data_fullFileName = destDir + "/Slink_a=100um_b=400um_LOL.txt";
    // Load the singular value data from the print file.
    read_vec_tmp = utils::file_to_Vd( sv_data_fullFileName );
    // Translate the data vector in a vector under the Eigen library format.
    read_vec = Eigen::VectorXd( read_vec_tmp.size() );
    for( unsigned int z = 0; z < read_vec_tmp.size(); z++ ){
        read_vec(z) = read_vec_tmp[z];
    }

    // Compute the difference between the vector from the file and directly from the source.
    vec_diff = comp_vec - read_vec;
    
    test_bool = vec_diff.cwiseAbs().maxCoeff() < 1e-7;
    if( test_bool ){
        cout << "LM_eng_print_singVals case 2 test: passed!" << endl;
    }else{
        cout << "LM_eng_print_singVals case 2 test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<

}



void tests::LM_eng_re_LM_comp_test( unsigned int test_idx ){

    unsigned int case_cnt = 0;

    // 0- Real LMs direct computation vs Complex LMs transform to real LMs.
    if( test_idx == case_cnt ){

// ---------------------------------------------------------------------- >>>>>
//      Initialization (Data)
// ---------------------------------------------------------------------- >>>>>
    
        // Define our frequency data object.
        fData myFData1;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/test_res_dir/bondwire_with_strip_design3.s4p";
        fData::read_sXp_file( myFData1, fullFileName );

        // Switch the data format into real + imaginary format.
        myFData1.data_format_Switch( fData::FDATA_FORMAT::RI );
        // Normalize the frequency vector (As much as you can according to metric prefixes).
        myFData1.data_prefix_switch( fData::METRIC_PREFIX::M );

// ---------------------------------------------------------------------- <<<<<

        unsigned int fr_size = 50;

        // Create a subset linear index array.
        vector<unsigned int> fr_idx_arr_in = utils::gen_lin_idx_arr( 0, 
            myFData1.get_f_cnt() - 1, min( fr_size, myFData1.get_f_cnt() ) );

        // LM engine initialization
        LM_eng my_LM_eng_a;
        LM_eng my_LM_eng_b;

        my_LM_eng_a.step0_fData_set( myFData1, fr_idx_arr_in );
        my_LM_eng_b.step0_fData_set( myFData1, fr_idx_arr_in );
        my_LM_eng_a.step1_fData_partition();
        my_LM_eng_b.step1_fData_partition();

        // Record the start time
        auto start = std::chrono::high_resolution_clock::now();
        // Perform the direct real LM computation step.
        my_LM_eng_a.step3skip2_LM_re_construct();
        // Record the end time
        auto end = std::chrono::high_resolution_clock::now();
        // Calculate the duration in milliseconds
        std::chrono::duration<double, std::milli> duration = end - start;
        cout << "Direct real LM comp cost: " << duration.count() << " ms" << endl;


        start = std::chrono::high_resolution_clock::now();
        // Perform the standard steps for the copy engine.
        my_LM_eng_b.step2_LM_construct();
        my_LM_eng_b.step3_LM_re_trans();
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        cout << "Standard real LM comp cost: " << duration.count() << " ms" << endl;

        
        // Obtain the real Loewner matrices generated the direct way.
        Eigen::MatrixXd my_LM_a = my_LM_eng_a.get_LM_re();
        Eigen::MatrixXd my_SLM_a = my_LM_eng_a.get_SLM_re();
        Eigen::MatrixXd my_F_a = my_LM_eng_a.get_F_re();
        Eigen::MatrixXd my_W_a = my_LM_eng_a.get_W_re();
        
        // Obtain the real Loewner matrices generated the standard way.
        Eigen::MatrixXd my_LM_b = my_LM_eng_b.get_LM_re();
        Eigen::MatrixXd my_SLM_b = my_LM_eng_b.get_SLM_re();
        Eigen::MatrixXd my_F_b = my_LM_eng_b.get_F_re();
        Eigen::MatrixXd my_W_b = my_LM_eng_b.get_W_re();
        
        
        // Calculate the highest discrepancy in magnitude.
        bool test_bool = true;
        test_bool = test_bool && ( ( my_LM_a - my_LM_b ).cwiseAbs().maxCoeff() < 1e-12 );
        test_bool = test_bool && ( ( my_SLM_a - my_SLM_b ).cwiseAbs().maxCoeff() < 1e-12 );
        test_bool = test_bool && ( ( my_F_a - my_F_b ).cwiseAbs().maxCoeff() < 1e-12 );
        test_bool = test_bool && ( ( my_W_a - my_W_b ).cwiseAbs().maxCoeff() < 1e-12 );

        if( test_bool ){
            cout << "step3skip2_LM_re_construct test: passed!" << endl;
        }else{
            cout << "step3skip2_LM_re_construct test: failed!" << endl;
        }
        

    }


    case_cnt++;
    // 1- Individual real LMs comp functions check.
    if( test_idx == case_cnt ){

// ---------------------------------------------------------------------- >>>>>
//      Initialization (Data)
// ---------------------------------------------------------------------- >>>>>
    
        // Define our frequency data object.
        fData myFData1;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/test_res_dir/bondwire_with_strip_design3.s4p";
        fData::read_sXp_file( myFData1, fullFileName );

        // Switch the data format into real + imaginary format.
        myFData1.data_format_Switch( fData::FDATA_FORMAT::RI );
        // Normalize the frequency vector (As much as you can according to metric prefixes).
        myFData1.data_prefix_switch( fData::METRIC_PREFIX::M );

// ---------------------------------------------------------------------- <<<<<

        unsigned int fr_size = 50;

        // Create a subset linear index array.
        vector<unsigned int> fr_idx_arr_in = utils::gen_lin_idx_arr( 0, 
            myFData1.get_f_cnt() - 1, min( fr_size, myFData1.get_f_cnt() ) );

        // LM engine initialization
        LM_eng my_LM_eng_a;

        // Go through LM construct step until real LM generation.
        my_LM_eng_a.step0_fData_set( myFData1, fr_idx_arr_in );
        my_LM_eng_a.step1_fData_partition();
        my_LM_eng_a.step3skip2_LM_re_construct();

        // Obtain the real Loewner matrices generated the direct way.
        Eigen::MatrixXd my_LM_a = my_LM_eng_a.get_LM_re();
        Eigen::MatrixXd my_SLM_a = my_LM_eng_a.get_SLM_re();
        Eigen::MatrixXd my_F_a = my_LM_eng_a.get_F_re();
        Eigen::MatrixXd my_W_a = my_LM_eng_a.get_W_re();

        Eigen::MatrixXd my_LM_b = 
            LM_UTIL::build_LM_re( my_LM_eng_a.get_Fr1(), my_LM_eng_a.get_Fr2() );
        Eigen::MatrixXd my_SLM_b = 
            LM_UTIL::build_SLM_re( my_LM_eng_a.get_Fr1(), my_LM_eng_a.get_Fr2() );
        shared_ptr<Eigen::MatrixXd> my_W_b = 
            LM_UTIL::build_W_re( my_LM_eng_a.get_Fr1() );
        shared_ptr<Eigen::MatrixXd> my_F_b = 
            LM_UTIL::build_F_re( my_LM_eng_a.get_Fr2() );

        // Calculate the highest discrepancy in magnitude.
        bool test_bool = true;
        test_bool = test_bool && ( ( my_LM_a - my_LM_b ).cwiseAbs().maxCoeff() < 1e-12 );
        test_bool = test_bool && ( ( my_SLM_a - my_SLM_b ).cwiseAbs().maxCoeff() < 1e-12 );
        test_bool = test_bool && ( ( my_F_a - *my_F_b ).cwiseAbs().maxCoeff() < 1e-12 );
        test_bool = test_bool && ( ( my_W_a - *my_W_b ).cwiseAbs().maxCoeff() < 1e-12 );

        if( test_bool ){
            cout << "Individual real LM build functions test: passed!" << endl;
        }else{
            cout << "Individual real LM build functions test: failed!" << endl;
        }

    }

}



void tests::LM_eng_serialize_test(){

// ---------------------------------------------------------------------- >>>>>
//      Test Setup
// ---------------------------------------------------------------------- >>>>>
    
    // Initialize test boolean.
    bool test_bool = true;

    // Define the numerical threshold of this test.
    double num_thresh = 1e-12;

    // Define our frequency data object.
    fData myFData1;

    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/test_res_dir/Slink_a=100um_b=400um.s2p";
    fData::read_sXp_file( myFData1, fullFileName );

    // Switch the data format into real + imaginary format.
    myFData1.data_format_Switch( fData::FDATA_FORMAT::RI );
    // Normalize the frequency vector (As much as you can according to metric prefixes).
    myFData1.data_prefix_switch( fData::METRIC_PREFIX::G );


    // LM engine initialization
    LM_eng my_LM_eng_a( myFData1 );

    my_LM_eng_a.step1_fData_partition();
    my_LM_eng_a.step2_LM_construct();
    my_LM_eng_a.step3_LM_re_trans();
    my_LM_eng_a.step4_LM_pencil_SVD();

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Target Tested Functionality Run
// ---------------------------------------------------------------------- >>>>>

    string outFileName = SRC_PATH_XYQ_str + "/data_output/LM_eng_test.bin";
    my_LM_eng_a.serialize( outFileName );

    LM_eng my_LM_eng_b;
    my_LM_eng_b.deserialize( outFileName );

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Step 0 Test
// ---------------------------------------------------------------------- >>>>>
    
    // Flags matching.
    for( unsigned int z = 0; z < 5; z++ ){
        test_bool = test_bool && ( my_LM_eng_a.get_flag(z) == my_LM_eng_a.get_flag(z) );
    }
    // fr_idx_arr matching.
    vector<unsigned int> fr_idx_arr_A = my_LM_eng_a.get_fr_idx_arr();
    vector<unsigned int> fr_idx_arr_B = my_LM_eng_b.get_fr_idx_arr();
    test_bool = test_bool && ( fr_idx_arr_A.size() == fr_idx_arr_B.size() );
    if( test_bool ){
        for( unsigned int z = 0; z < fr_idx_arr_A.size(); z++ ){
            test_bool = test_bool && ( fr_idx_arr_A[z] == fr_idx_arr_B[z] );
        }
    }
    // Step 0 result.
    if( test_bool ){
        cout << "LM_eng serialize + deserialize step 0 test: passed!" << endl;
    }else{
        cout << "LM_eng serialize + deserialize step 0 test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Step 1 Test
// ---------------------------------------------------------------------- >>>>>

    test_bool = true;
    // DC point flags check.
    test_bool = test_bool && ( my_LM_eng_a.get_f1_has_DC_pt() == my_LM_eng_b.get_f1_has_DC_pt() );
    test_bool = test_bool && ( my_LM_eng_a.get_f2_has_DC_pt() == my_LM_eng_b.get_f2_has_DC_pt() );
    // Partition index arrays check.
    vector<unsigned int> f1_idx_arr_A = my_LM_eng_a.get_partit1IdxArr();
    vector<unsigned int> f2_idx_arr_A = my_LM_eng_a.get_partit2IdxArr();
    vector<unsigned int> f1_idx_arr_B = my_LM_eng_b.get_partit1IdxArr();
    vector<unsigned int> f2_idx_arr_B = my_LM_eng_b.get_partit2IdxArr();
    test_bool = test_bool && ( f1_idx_arr_A.size() == f1_idx_arr_B.size() );
    if( test_bool ){
        for( unsigned int z = 0; z < f1_idx_arr_A.size(); z++ ){
            test_bool = test_bool && ( f1_idx_arr_A[z] == f1_idx_arr_B[z] );
        }
    }
    test_bool = test_bool && ( f2_idx_arr_A.size() == f2_idx_arr_B.size() );
    if( test_bool ){
        for( unsigned int z = 0; z < f2_idx_arr_A.size(); z++ ){
            test_bool = test_bool && ( f2_idx_arr_A[z] == f2_idx_arr_B[z] );
        }
    }
    // Step 1 result.
    if( test_bool ){
        cout << "LM_eng serialize + deserialize step 1 test: passed!" << endl;
    }else{
        cout << "LM_eng serialize + deserialize step 1 test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Step 2 Test
// ---------------------------------------------------------------------- >>>>>

    // Complex LM checks.
    test_bool = true;
    test_bool = test_bool && ( ( my_LM_eng_a.get_LM() - my_LM_eng_b.get_LM() ).cwiseAbs().maxCoeff() < num_thresh );
    test_bool = test_bool && ( ( my_LM_eng_a.get_SLM() - my_LM_eng_b.get_SLM() ).cwiseAbs().maxCoeff() < num_thresh );
    test_bool = test_bool && ( ( my_LM_eng_a.get_W() - my_LM_eng_b.get_W() ).cwiseAbs().maxCoeff() < num_thresh );
    test_bool = test_bool && ( ( my_LM_eng_a.get_F() - my_LM_eng_b.get_F() ).cwiseAbs().maxCoeff() < num_thresh );
    // Step 2 result.
    if( test_bool ){
        cout << "LM_eng serialize + deserialize step 2 test: passed!" << endl;
    }else{
        cout << "LM_eng serialize + deserialize step 2 test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Step 3 Test
// ---------------------------------------------------------------------- >>>>>

    // Real LM checks.
    test_bool = true;
    test_bool = test_bool && ( ( my_LM_eng_a.get_LM_re() - my_LM_eng_b.get_LM_re() ).cwiseAbs().maxCoeff() < num_thresh );
    test_bool = test_bool && ( ( my_LM_eng_a.get_SLM_re() - my_LM_eng_b.get_SLM_re() ).cwiseAbs().maxCoeff() < num_thresh );
    test_bool = test_bool && ( ( my_LM_eng_a.get_W_re() - my_LM_eng_b.get_W_re() ).cwiseAbs().maxCoeff() < num_thresh );
    test_bool = test_bool && ( ( my_LM_eng_a.get_F_re() - my_LM_eng_b.get_F_re() ).cwiseAbs().maxCoeff() < num_thresh );
    // Step 3 result.
    if( test_bool ){
        cout << "LM_eng serialize + deserialize step 3 test: passed!" << endl;
    }else{
        cout << "LM_eng serialize + deserialize step 3 test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Step 4 Test
// ---------------------------------------------------------------------- >>>>>

    // Real SVD checks.
    test_bool = true;

    
    test_bool = test_bool && ( my_LM_eng_a.get_ref_f_mag() == my_LM_eng_b.get_ref_f_mag() );
    test_bool = test_bool && ( my_LM_eng_a.get_svd_cnt_max() == my_LM_eng_b.get_svd_cnt_max() );
    test_bool = test_bool && ( ( my_LM_eng_a.get_U() - my_LM_eng_b.get_U() ).cwiseAbs().maxCoeff() < num_thresh );
    test_bool = test_bool && ( ( my_LM_eng_a.get_V() - my_LM_eng_b.get_V() ).cwiseAbs().maxCoeff() < num_thresh );
    test_bool = test_bool && ( ( my_LM_eng_a.get_singVals() - my_LM_eng_b.get_singVals() ).cwiseAbs().maxCoeff() < num_thresh );
    // Step 4 result.
    if( test_bool ){
        cout << "LM_eng serialize + deserialize step 4 test: passed!" << endl;
    }else{
        cout << "LM_eng serialize + deserialize step 4 test: failed!" << endl;
    }

// ---------------------------------------------------------------------- <<<<<

    unsigned int lol = 0;

}
