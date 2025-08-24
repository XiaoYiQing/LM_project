#include "LM_eng.h"



const double LM_eng::NUM_THRESH = 1.0e-12;


// ====================================================================== >>>>>
//      Static Support Functions
// ====================================================================== >>>>>

shared_ptr<LM_eng> LM_eng::print_singVals( const string& fullFileName, 
    const string& destDir ){

    // Obtain the specific parts of the full file name (dir, stem, ext).
    std::filesystem::path fullFilePath( fullFileName );
    std::string fileDir = fullFilePath.parent_path().string();
    std::string fileStem = fullFilePath.stem().string();
    std::string fileExt = fullFilePath.extension().string();
    // Abort if any element is missing, abort the operation.
    if( fileDir.size() == 0 || fileStem.size() == 0 || fileExt.size() == 0 ){
        throw std::invalid_argument( "Full file name incomplete." );
    }

    /*
    Regex for determining the positive integer value X in the pattern ".sXp".
    */
    regex pattern(R"(\.s(\d+)p)");
    // The match result variable.
    smatch matches;

    // Define our frequency data object.
    fData myFData;

    // Parse LTspice text data output format.
    if( fileExt == ".txt" ){
        
        try{
            fData::read_LTspice_Sp_file( myFData, fullFileName );
        }catch( const std::invalid_argument& e ){
            std::cerr << "print_singVals aborted: " << e.what() << '\n';
            shared_ptr<LM_eng> tmp;
            return tmp; 
        }catch( const std::runtime_error& e ){
            std::cerr << "print_singVals aborted: " << e.what() << '\n';
            shared_ptr<LM_eng> tmp;
            return tmp; 
        }

    // Parse touchstone files.
    }else if( regex_match( fileExt, matches, pattern ) ){

        // Try to obtain data from specified data file.
        try{
            // Obtain the data from the target data file and insert into the fData object.
            fData::read_sXp_file( myFData, fullFileName );
        }catch( const std::invalid_argument& e ){
            std::cerr << "print_singVals aborted: " << e.what() << '\n';
            shared_ptr<LM_eng> tmp;
            return tmp; 
        }catch( const std::runtime_error& e ){
            std::cerr << "print_singVals aborted: " << e.what() << '\n';
            shared_ptr<LM_eng> tmp;
            return tmp; 
        }

    // Unrecognized extension case.
    }else{
        throw std::invalid_argument( "print_singVals aborted: target data file has unrecognized extension." );
    }

    
    
    // Switch the data format into real + imaginary format.
    myFData.data_format_Switch( fData::FDATA_FORMAT::RI );
    // Normalize the frequency vector (As much as you can according to metric prefixes).
    myFData.data_prefix_switch( fData::METRIC_PREFIX::M );

    
    // LM engine initialization.
    shared_ptr<LM_eng> myEng = make_shared<LM_eng>( myFData );
    try{
        myEng->step1_fData_partition();
        myEng->step2_LM_construct();
        myEng->step3_LM_re_trans();
        myEng->step4_LM_pencil_SVD();
    }catch( const std::runtime_error& e ){
        std::cerr << "print_singVals aborted: " << e.what() << '\n';
        shared_ptr<LM_eng> tmp;
        return tmp; 
    }

    string dataFileStem = fileStem + "_sv";

    utils::vec_to_file( destDir, dataFileStem, myEng->get_singVals(), 0 );

    return myEng;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Constructor
// ====================================================================== >>>>>

LM_eng::LM_eng(){

}

LM_eng::LM_eng( const fData& inData ){

    this->myFData = inData;
    flag0_data_set = true;

}

// ====================================================================== <<<<<




// ====================================================================== >>>>>
//      Major LM System Steps
// ====================================================================== >>>>>

void LM_eng::step1_fData_partition(){

    if( !flag0_data_set ){
        throw::runtime_error( "Step 1 cannot be executed: step 0 not set (starting data insertion)." );
    }

    // Create a subset linear index array.
    vector< unsigned int > fr_idx_arr = 
        utils::gen_lin_idx_arr( 0, this->myFData.get_f_cnt() - 1, s1_fr_len );

    // Reset the reduced frequency set variable.
    this->myFr.reset();
    // Create a fData subset.
    this->myFr = this->myFData.red_partit( fr_idx_arr );

    // Generate the partition index arrays.
    vector< vector< unsigned int > > index_arrs = this->myFr->gen_2_partit_idx_arr();
    this->partit1IdxArr = index_arrs.at(0);
    this->partit2IdxArr = index_arrs.at(1);

    // Check for DC point in either partitions.
    this->f1_has_DC_pt = myFr->get_fval_at( partit1IdxArr.at(0) ) == 0;
    this->f2_has_DC_pt = myFr->get_fval_at( partit2IdxArr.at(0) ) == 0;

    // Set the tracking flag for step 1.
    this->flag1_data_prep = true;

}

void LM_eng::step2_LM_construct(){

    if( !flag1_data_prep ){
        throw::runtime_error( "Step 2 cannot be executed: step 1 not set (data prep)." );
    }

    // Create a fData partitions.
    shared_ptr<fData> myFr1 = this->myFr->red_partit( this->partit1IdxArr );
    shared_ptr<fData> myFr2 = this->myFr->red_partit( this->partit2IdxArr );

    // Generate the two partitions with their complex conjugates inserted 
    // in interleaving fashion.
    shared_ptr<fData> myFrc1 = myFr1->gen_cplx_conj_comb();
    shared_ptr<fData> myFrc2 = myFr2->gen_cplx_conj_comb();

    // Construct the Loewner Matrix using the two cconj injected partitions.
    this->LM = *LM_UTIL::build_LM( *myFrc1, *myFrc2 );
    // Construct the Loewner Matrix using the two cconj injected partitions.
    this->SLM = *LM_UTIL::build_SLM( *myFrc1, *myFrc2 );
    // Construct the W matrix vector using partition 1.
    this->W = *LM_UTIL::build_W( *myFrc1 );
    // Construct the F matrix vector using partition 2.
    this->F = *LM_UTIL::build_F( *myFrc2 );

    // Set the tracking flag for step 1.
    this->flag2_LM_const = true;

}


void LM_eng::step3_LM_re_trans(){
    
    if( !flag2_LM_const ){
        throw::runtime_error( "Step 3 cannot be executed: step 2 not set (LM construction)." );
    }

    // Partition sizes before cconj injection.
    unsigned int out_cnt = this->myFr->get_out_cnt();
    unsigned int in_cnt = this->myFr->get_in_cnt();
    unsigned int fr1_len = this->partit1IdxArr.size();
    unsigned int fr2_len = this->partit2IdxArr.size();

    // Build the left and right transformation matrices.
    Eigen::MatrixXcd myTMat_L = *LM_UTIL::build_reT_mat( this->f2_has_DC_pt, out_cnt, fr1_len );
    Eigen::MatrixXcd myTMat_R = *LM_UTIL::build_reT_mat( this->f1_has_DC_pt, out_cnt, fr2_len );
    // Obtain the hermitian of the right transform matrix.
    Eigen::MatrixXcd myTMat_L_herm = myTMat_L.conjugate().transpose();

    // Perform the transformation.
    Eigen::MatrixXcd myLM_re_tmp = ( myTMat_L_herm*this->LM )*myTMat_R;
    Eigen::MatrixXcd mySLM_re_tmp = ( myTMat_L_herm*this->SLM )*myTMat_R;
    Eigen::MatrixXcd myW_re_tmp = this->W*myTMat_R;
    Eigen::MatrixXcd myF_re_tmp = myTMat_L_herm*this->F;

    // Check for real matrices.
    bool match_bool = true;
    match_bool = match_bool && ( myLM_re_tmp.imag().cwiseAbs().maxCoeff() < LM_eng::NUM_THRESH );
    match_bool = match_bool && ( mySLM_re_tmp.imag().cwiseAbs().maxCoeff() < LM_eng::NUM_THRESH );
    match_bool = match_bool && ( myW_re_tmp.imag().cwiseAbs().maxCoeff() < LM_eng::NUM_THRESH );
    match_bool = match_bool && ( myF_re_tmp.imag().cwiseAbs().maxCoeff() < LM_eng::NUM_THRESH );
    if( !match_bool ){
        cout << "WARNING: SFLM real matrices check has failed. Non-negligeable imaginary part remnant detected." << endl;
    }

    // Obtain purely real defintion of the matrices.
    this->LM_re = myLM_re_tmp.real();
    this->SLM_re = mySLM_re_tmp.real();
    this->W_re = myW_re_tmp.real();
    this->F_re = myF_re_tmp.real();

    // Generate a random test point and evaluate the full LM transfer function.
    // if( !f2_has_DC_pt && !f1_has_DC_pt ){
    //     unsigned int test_f_idx = utils::rIntGen( 0, this->myFData.get_f_cnt() - 1, 1 )->at(0);
    //     complex<double> test_f = this->myFData.get_cplx_f_at( test_f_idx );
    //     Eigen::MatrixXcd tmpAns = 
    //         this->W_re*( ( - test_f*this->LM_re + this->SLM_re ).inverse() )*this->F_re;
    //     Eigen::MatrixXcd ansDiff = this->myFData.get_cplxData_at_f( test_f_idx ) - tmpAns;
    //     match_bool = true;
    //     match_bool = match_bool && ( ansDiff.cwiseAbs2().maxCoeff() < 1e-12 );
    //     cout << "Full sized LM system evaluation test (Not mandatory to pass): " << match_bool << endl;
    // }
    
    this->flag3_re_trans = true;

}


void LM_eng::step4_LM_pencil_SVD(){

    if( !flag3_re_trans ){
        throw::runtime_error( "Step 4 cannot be executed: step 3 not set (LM real transform)." );
    }

    // Obtain a reference frequency value.
    this->ref_f_mag = this->myFr->get_fval_at( this->myFr->get_f_cnt() - 1 );
    
    // Construct the LM pencil.
    shared_ptr<Eigen::MatrixXd> LM_pen = 
        LM_UTIL::build_LM_pencil( this->ref_f_mag, this->LM_re, this->SLM_re );

    // Perform SVD.
    Eigen::JacobiSVD<Eigen::MatrixXd> svdResObj( *LM_pen, Eigen::ComputeFullU | Eigen::ComputeFullV );
    // Get the singular values
    this->singVals = svdResObj.singularValues();
    // Get the left singular vectors (U)
    this->U = svdResObj.matrixU();
    // Get the right singular vectors (V)
    this->V = svdResObj.matrixV();

    flag4_pen_SVD = true;

}


shared_ptr<LTI_descSyst> LM_eng::step5_LM_to_tf( unsigned int svd_ret_cnt ){

    if( !flag4_pen_SVD ){
        throw::runtime_error( "Step 5 cannot be executed: step 4 not set (LM pencil SVD)." );
    }

    if( svd_ret_cnt > (unsigned int) singVals.size() ){
        throw::out_of_range( "Specified singular value index is out of range of available singular values." );
    }

    // Define number of outputs.
    unsigned int out_cnt = this->myFr->get_out_cnt();

    // Eigen::VectorXd singVals_r = this->singVals.segment( 0, svd_ret_cnt );
    
    Eigen::MatrixXd U_r = this->U.block( 0, 0, U.rows(), svd_ret_cnt );
    Eigen::MatrixXd V_r = this->V.block( 0, 0, V.rows(), svd_ret_cnt );

    // Perform the model reduction to obtain usable E, A, B, C matrices.
    Eigen::MatrixXd E_n = -1*( U_r.transpose() * this->LM_re * V_r );
    Eigen::MatrixXd A_n = -1*( U_r.transpose() * this->SLM_re * V_r );
    Eigen::MatrixXd C_n = this->W_re * V_r;
    Eigen::MatrixXd B_n = U_r.transpose() * this->F_re;
    Eigen::MatrixXd D_n = Eigen::MatrixXd::Zero( out_cnt, out_cnt );

    // Model generation.
    shared_ptr<LTI_descSyst> mySyst = 
        make_shared<LTI_descSyst>( E_n, A_n, B_n, C_n, D_n );

    return mySyst;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Function
// ====================================================================== >>>>>

bool LM_eng::get_flag( unsigned int flagIdx ) const{

    switch( flagIdx ){
    case 0:
        return this->flag0_data_set;
    case 1:
        return this->flag1_data_prep;
    case 2:
        return this->flag2_LM_const;
    case 3:
        return this->flag3_re_trans;
    case 4:
        return this->flag4_pen_SVD;
    default:
        throw std::invalid_argument( "Specified flag index does not exist." );
    };

}

void LM_eng::set_fData( const fData& inData ){

    // Set the engine's data with the given data.
    this->myFData = inData;

    // Reset the flags.
    this->flag0_data_set = false;
    this->flag1_data_prep = false;
    this->flag2_LM_const = false;
    this->flag3_re_trans = false;
    this->flag4_pen_SVD = false;

}

double LM_eng::get_ref_f_mag() const{

    if( !this->flag4_pen_SVD ){
        throw std::runtime_error( "Cannot return LM: step4 (LM pencil SVD) has not been set." );
    }
    return this->ref_f_mag;

}


unsigned int LM_eng::get_out_cnt() const{
    if( !this->flag1_data_prep ){
        throw std::runtime_error( "Cannot return output count: step1 (data preparation) has not been set." );
    }
    return this->myFData.get_out_cnt();
}
unsigned int LM_eng::get_in_cnt() const{
    if( !this->flag1_data_prep ){
        throw std::runtime_error( "Cannot return output count: step1 (data preparation) has not been set." );
    }
    return this->myFData.get_in_cnt();
}



shared_ptr<fData> LM_eng::get_Fr1() const{
    if( !this->flag1_data_prep ){
        throw std::runtime_error( "Cannot return reduced frequency partition 1: step1 (data preparation) has not been set." );
    }
    return this->myFr->red_partit( this->partit1IdxArr );
}
shared_ptr<fData> LM_eng::get_Fr2() const{
    if( !this->flag1_data_prep ){
        throw std::runtime_error( "Cannot return reduced frequency partition 2: step1 (data preparation) has not been set." );
    }
    return this->myFr->red_partit( this->partit2IdxArr );
}

shared_ptr<fData> LM_eng::get_Frc1() const{
    if( !this->flag1_data_prep ){
        throw std::runtime_error( "Cannot return reduced complex conjugate frequency partition 1: step1 (data preparation) has not been set." );
    }
    return this->myFr->red_partit( this->partit1IdxArr )->gen_cplx_conj_comb();
}
shared_ptr<fData> LM_eng::get_Frc2() const{
    if( !this->flag1_data_prep ){
        throw std::runtime_error( "Cannot return reduced complex conjugate frequency partition 2: step1 (data preparation) has not been set." );
    }
    return this->myFr->red_partit( this->partit2IdxArr )->gen_cplx_conj_comb();
}


Eigen::MatrixXcd LM_eng::get_LM() const{
    if( !this->flag2_LM_const ){
        throw std::runtime_error( "Cannot return LM: step2 (LM construction) has not been set." );
    }
    return this->LM;
}
Eigen::MatrixXcd LM_eng::get_SLM() const{
    if( !this->flag2_LM_const ){
        throw std::runtime_error( "Cannot return SLM: step2 (LM construction) has not been set." );
    }
    return this->SLM;
}
Eigen::MatrixXcd LM_eng::get_W() const{
    if( !this->flag2_LM_const ){
        throw std::runtime_error( "Cannot return W: step2 (LM construction) has not been set." );
    }
    return this->W;
}
Eigen::MatrixXcd LM_eng::get_F() const{
    if( !this->flag2_LM_const ){
        throw std::runtime_error( "Cannot return F: step2 (LM construction) has not been set." );
    }
    return this->F;
}


Eigen::MatrixXcd LM_eng::get_LM_re() const{
    if( !this->flag3_re_trans){
        throw std::runtime_error( "Cannot return real LM: step3 (LM real transform) has not been set." );
    }
    return this->LM_re;
}
Eigen::MatrixXcd LM_eng::get_SLM_re() const{
    if( !this->flag3_re_trans ){
        throw std::runtime_error( "Cannot return real SLM: step3 (LM real transform) has not been set." );
    }
    return this->SLM_re;
}
Eigen::MatrixXcd LM_eng::get_W_re() const{
    if( !this->flag3_re_trans ){
        throw std::runtime_error( "Cannot return real W: step3 (LM real transform) has not been set." );
    }
    return this->W_re;
}
Eigen::MatrixXcd LM_eng::get_F_re() const{
    if( !this->flag3_re_trans ){
        throw std::runtime_error( "Cannot return real F: step3 (LM real transform) has not been set." );
    }
    return this->F_re;
}



Eigen::VectorXd LM_eng::get_singVals() const{

    if( !this->flag4_pen_SVD ){
        throw std::runtime_error( "Cannot return singular values: step4 (LM pencil SVD) has not been set." );
    }
    return this->singVals;

}

// Obtain the left singular vectors generated from the current LM pencil.
Eigen::MatrixXd LM_eng::get_U() const{
    if( !this->flag4_pen_SVD ){
        throw std::runtime_error( "Cannot return left singular vectors: step4 (LM pencil SVD) has not been set." );
    }
    return this->U;
}
// Obtain the right singular vectors generated from the current LM pencil.
Eigen::MatrixXd LM_eng::get_V() const{
    if( !this->flag4_pen_SVD ){
        throw std::runtime_error( "Cannot return right singular vectors: step4 (LM pencil SVD) has not been set." );
    }
    return this->V;
}

// ====================================================================== <<<<<



shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_LM( const fData& f1Data, const fData& f2Data ){

    // Obtain the size of the two partitions.
    unsigned int f1Size = f1Data.get_f_cnt();
    unsigned int f2Size = f2Data.get_f_cnt();
    if( f1Size == 0 || f2Size == 0 ){
        throw std::invalid_argument( "Empty frequency data inputs are not allowed for constructing the Loewner Matrix." );
    }
    
    // Obtain the number of outputs and inputs.
    unsigned int out_cnt = f1Data.get_out_cnt();
    unsigned int in_cnt = f1Data.get_in_cnt();

    // Determine the expected size of the Loewner Matrix.
    unsigned int row_cnt = f2Size*out_cnt;
    unsigned int col_cnt = f1Size*in_cnt;

    // Initialize the Loewner Matrix.
    shared_ptr<Eigen::MatrixXcd> LM = std::make_shared<Eigen::MatrixXcd>( row_cnt, col_cnt );

    // Initialize temporary matrix representing the current sub-block being computed.
    Eigen::MatrixXcd LM_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Initialize current freq. data from each partition.
    Eigen::MatrixXcd f1_D_j = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd f2_D_i = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Initialize current frequencies from each partition.
    complex<double> f1_j = complex<double>(0,0); 
    complex<double> f2_i = complex<double>(0,0);
    

    // Initialize current lead coordinate where the sub-block is to be inserted in the
    // Loewner Matrix.
    unsigned int lead_x = 0, lead_y = 0;

    // Build the Loewner Matrix block by block.
    for( unsigned int i = 0; i < f2Size; i++ ){

        // Obtain the current f2 freq. and data.
        f2_i = f2Data.get_cplx_f_at( i );
        f2_D_i = f2Data.get_cplxData_at_f( i );

        for( unsigned int j = 0; j < f1Size; j++ ){

            // Obtain the current f1 freq. and data.
            f1_j = f1Data.get_cplx_f_at( j );
            f1_D_j = f1Data.get_cplxData_at_f( j );

            // Obtain the matrix coordinate of the upper-left leading point of the current
            // LM sub-block being computed.
            lead_x = i*out_cnt;
            lead_y = j*in_cnt;

            // Compute the current Loewner Matrix block.
            LM_ij = ( f2_D_i - f1_D_j )/( f2_i - f1_j );

            // Insert the current calculated block into its part in the full matrix.
            LM->block( lead_x, lead_y, out_cnt, in_cnt) = LM_ij;

        }

    }


    return LM;

}



shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_SLM( const fData& f1Data, const fData& f2Data ){

    // Obtain the size of the two partitions.
    unsigned int f1Size = f1Data.get_f_cnt();
    unsigned int f2Size = f2Data.get_f_cnt();
    if( f1Size == 0 || f2Size == 0 ){
        throw std::invalid_argument( "Empty frequency data inputs are not allowed for constructing the shifted Loewner Matrix." );
    }

    // Obtain the number of outputs and inputs.
    unsigned int out_cnt = f1Data.get_out_cnt();
    unsigned int in_cnt = f1Data.get_in_cnt();

    // Determine the expected size of the shifted Loewner Matrix.
    unsigned int row_cnt = f2Size*out_cnt;
    unsigned int col_cnt = f1Size*in_cnt;

    // Initialize the shifted Loewner Matrix.
    shared_ptr<Eigen::MatrixXcd> sLM = std::make_shared<Eigen::MatrixXcd>( row_cnt, col_cnt );

    // Initialize temporary matrix representing the current sub-block being computed.
    Eigen::MatrixXcd sLM_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Initialize current freq. data from each partition.
    Eigen::MatrixXcd f1_D_j = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd f2_D_i = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Initialize current frequencies from each partition.
    complex<double> f1_j = complex<double>(0,0); 
    complex<double> f2_i = complex<double>(0,0);

    // Initialize current lead coordinate where the sub-block is to be inserted in the
    // shifted Loewner Matrix.
    unsigned int lead_x = 0, lead_y = 0;

    // Build the shifted Loewner Matrix block by block.
    for( unsigned int i = 0; i < f2Size; i++ ){

        // Obtain the current f2 freq. and data.
        f2_i = f2Data.get_cplx_f_at( i );
        f2_D_i = f2Data.get_cplxData_at_f( i );

        for( unsigned int j = 0; j < f1Size; j++ ){

            // Obtain the current f1 freq. and data.
            f1_j = f1Data.get_cplx_f_at( j );
            f1_D_j = f1Data.get_cplxData_at_f( j );

            // Obtain the matrix coordinate of the upper-left leading point of the current
            // sLM sub-block being computed.
            lead_x = i*out_cnt;
            lead_y = j*in_cnt;

            // Compute the current shifted Loewner Matrix block.
            sLM_ij = ( f2_i*f2_D_i - f1_j*f1_D_j )/( f2_i - f1_j );

            // Insert the current calculated block into its part in the full matrix.
            sLM->block( lead_x, lead_y, out_cnt, in_cnt) = sLM_ij;

        }

    }

    return sLM;

}

shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_W( const fData& f1Data ){

    // Obtain the size of the two partitions.
    unsigned int f1Size = f1Data.get_f_cnt();
    if( f1Size == 0 ){
        throw std::invalid_argument( "Empty frequency data input is not allowed for constructing the W matrix vector." );
    }

    // Obtain the number of outputs and inputs.
    unsigned int out_cnt = f1Data.get_out_cnt();
    unsigned int in_cnt = f1Data.get_in_cnt();

    // Determine the expected size of the shifted Loewner Matrix.
    unsigned int col_cnt = f1Size*in_cnt;

    // Initialize the W matrix vector.
    shared_ptr<Eigen::MatrixXcd> W = std::make_shared<Eigen::MatrixXcd>( out_cnt, col_cnt );

    // Initialize current freq. data from each partition.
    Eigen::MatrixXcd f1_D_j = Eigen::MatrixXcd( out_cnt, in_cnt );

    // Initialize current lead coordinate where the sub-block is to be inserted in the
    // shifted Loewner Matrix.
    unsigned int lead_x = 0, lead_y = 0;

    for( unsigned int j = 0; j < f1Size; j++ ){

        // Obtain the current f1 freq. and data.
        f1_D_j = f1Data.get_cplxData_at_f( j );

        // Obtain the matrix coordinate of the upper-left leading point of the current
        // sLM sub-block being computed.
        lead_y = j*in_cnt;

        // Insert the current calculated block into its part in the full matrix.
        W->block( lead_x, lead_y, out_cnt, in_cnt) = f1_D_j;

    }

    return W;

}


shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_F( const fData& f2Data ){

    // Obtain the size of the two partitions.
    unsigned int f2Size = f2Data.get_f_cnt();
    if( f2Size == 0 ){
        throw std::invalid_argument( "Empty frequency data input is not allowed for constructing the F matrix vector." );
    }

    // Obtain the number of outputs and inputs.
    unsigned int out_cnt = f2Data.get_out_cnt();
    unsigned int in_cnt = f2Data.get_in_cnt();

    // Determine the expected size of the F matrix vector.
    unsigned int row_cnt = f2Size*out_cnt;

    // Initialize the W matrix vector.
    shared_ptr<Eigen::MatrixXcd> F = std::make_shared<Eigen::MatrixXcd>( row_cnt, in_cnt );

    // Initialize current freq. data from each partition.
    Eigen::MatrixXcd f2_D_i = Eigen::MatrixXcd( out_cnt, in_cnt );

    // Initialize current lead coordinate where the sub-block is to be inserted in the
    // shifted Loewner Matrix.
    unsigned int lead_x = 0, lead_y = 0;

    for( unsigned int i = 0; i < f2Size; i++ ){

        // Obtain the current f1 freq. and data.
        f2_D_i = f2Data.get_cplxData_at_f( i );

        // Obtain the matrix coordinate of the upper-left leading point of the current
        // sLM sub-block being computed.
        lead_x = i*out_cnt;

        // Insert the current calculated block into its part in the full matrix.
        F->block( lead_x, lead_y, out_cnt, in_cnt) = f2_D_i;

    }

    return F;

}


shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_LM_pencil( complex<double> ref_f, const Eigen::MatrixXcd& LM, const Eigen::MatrixXcd& SLM ){

    unsigned int row_cnt = LM.rows();
    unsigned int col_cnt = LM.cols();
    if( row_cnt != SLM.rows() || col_cnt != SLM.cols() ){
        throw std::invalid_argument( "The LM and the SLM must shared the same dimensions." );
    }
    
    shared_ptr<Eigen::MatrixXcd> LM_pen = std::make_shared<Eigen::MatrixXcd>( row_cnt, col_cnt );

    *LM_pen = ref_f*LM - SLM;

    return LM_pen;

}


shared_ptr<Eigen::MatrixXd> LM_UTIL::build_LM_pencil( double ref_f, const Eigen::MatrixXd& LM, 
    const Eigen::MatrixXd& SLM )
{
    unsigned int row_cnt = LM.rows();
    unsigned int col_cnt = LM.cols();
    if( row_cnt != SLM.rows() || col_cnt != SLM.cols() ){
        throw std::invalid_argument( "The LM and the SLM must shared the same dimensions." );
    }

    shared_ptr<Eigen::MatrixXd> LM_pen = std::make_shared<Eigen::MatrixXd>( row_cnt, col_cnt );

    *LM_pen = ref_f*LM - SLM;

    return LM_pen;

}


shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_reT_mat( bool has_DC_pt, unsigned int sub_mat_size, unsigned int sub_blk_cnt ){


    // Obtain the total number of rows and columns of the final transformation matrix.
    unsigned int row_cnt = 2*sub_blk_cnt*sub_mat_size;
    unsigned int col_cnt = 2*sub_blk_cnt*sub_mat_size;

    if( has_DC_pt ){
        row_cnt -= sub_mat_size;
        col_cnt -= sub_mat_size;    
    }

    // Complex 1.
    complex<double> im1(0,1);
    // Initialize the standard identity matrix.
    Eigen::MatrixXd I_mat = Eigen::MatrixXd::Identity(sub_mat_size, sub_mat_size);
    // Initialize the unit sub-block transformation matrix.
    Eigen::MatrixXcd T_unit = Eigen::MatrixXcd( 2*sub_mat_size, 2*sub_mat_size );
    T_unit.block( 0, 0, sub_mat_size, sub_mat_size ) = I_mat;
    T_unit.block( 0, sub_mat_size, sub_mat_size, sub_mat_size ) = -im1*I_mat;
    T_unit.block( sub_mat_size, 0, sub_mat_size, sub_mat_size ) = I_mat;
    T_unit.block( sub_mat_size, sub_mat_size, sub_mat_size, sub_mat_size ) = im1*I_mat;
    T_unit *= 1/( std::sqrt(2.0) );

    // Initialize the transformation matrix.
    shared_ptr<Eigen::MatrixXcd> T_mat = std::make_shared<Eigen::MatrixXcd>( row_cnt, col_cnt );
    T_mat->setZero();
    
    // Initialize current sub-block's lead coordinate.
    unsigned int lead_x = 0;
    // Initialize iteration index.
    unsigned int z = 0;
    // In case of DC point, fill the first block, then increment the iteration index.
    if( has_DC_pt ){
        T_mat->block( lead_x, lead_x, sub_mat_size, sub_mat_size ) = I_mat;
        lead_x += sub_mat_size;
        z++;
    }

    for( ; z < sub_blk_cnt; z++ ){

        // Insert the current sub-block.
        T_mat->block( lead_x, lead_x, 2*sub_mat_size, 2*sub_mat_size ) = T_unit;

        // Increment leading point on the matrix diagonal.
        lead_x += 2*sub_mat_size;

    }

    return T_mat;

}

