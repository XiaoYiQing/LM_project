#include "LM_eng.h"



const double LM_eng::NUM_THRESH = 1.0e-12;

const unsigned int LM_eng::STD_RED_FSET_SIZE = 100;

// ====================================================================== >>>>>
//      Data Printing Functions
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
    myFData.f_rescale();

    
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


void LM_eng::print_singVals( const LM_eng& tar_LM_eng, const string& fileStem, 
    const string& destDir )
{

    try{
        utils::vec_to_file( destDir, fileStem, tar_LM_eng.get_singVals(), 0 );
    }catch( const runtime_error& e ){
        cerr << "print_singVals aborted: " << e.what() << endl;
    }catch( const exception& e ){
        cerr << "print_singVals aborted due to unexpected error: " << e.what() << endl;
    } 

    cout << "print_singVals-> succesfully written singular values to: ";
    cout << destDir << "/" << fileStem << ".txt" << endl;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Constructor
// ====================================================================== >>>>>

LM_eng::LM_eng(){

}

LM_eng::LM_eng( const fData& inData ){

    this->step0_fData_set( inData );

}

// ====================================================================== <<<<<




// ====================================================================== >>>>>
//      Major LM System Steps
// ====================================================================== >>>>>

void LM_eng::step0_fData_set( const fData& inData ){

    // Create a subset linear index array.
    vector<unsigned int> fr_idx_arr_in = 
        utils::gen_lin_idx_arr( 0, inData.get_f_cnt() - 1, min( STD_RED_FSET_SIZE, inData.get_f_cnt() ) );
    
    // Follow the standard step 0 procedure.
    step0_fData_set( inData, fr_idx_arr_in );

}

void LM_eng::step0_fData_set( const fData& inData, const vector<unsigned int>& fr_idx_arr_in ){

    // Set the engine's data with the given data.
    this->myFData = inData;

    vector<unsigned int> fr_idx_arr_in_tmp = fr_idx_arr_in;
    // Sort the input index vector in ascending order.
    std::sort( fr_idx_arr_in_tmp.begin(), fr_idx_arr_in_tmp.end() );

    // Check for repeated indices.
    for( unsigned int z = 0; z < fr_idx_arr_in_tmp.size() - 1; z++ ){
        if( fr_idx_arr_in_tmp.at(z) == fr_idx_arr_in_tmp.at(z+1) ){
            throw std::invalid_argument( "Repeated entries in the reduced frequency set index vector." );
        }
    }
    // Check for out of bound indexing.
    if( fr_idx_arr_in_tmp.at( fr_idx_arr_in_tmp.size()-1 >= myFData.get_f_cnt() ) ){
        throw std::out_of_range( "Reduced frequency set index vector has indices exceeding number of available frequency entries." );
    }

    // Assign the local reduced set index vector.
    this->fr_idx_arr = fr_idx_arr_in_tmp;

    // Reset the flags.
    this->flag0_data_set = true;
    this->flag1_data_prep = false;
    this->flag2_LM_const = false;
    this->flag3_re_trans = false;
    this->flag4_pen_SVD = false;

}

void LM_eng::step1_fData_partition(){

    if( !flag0_data_set ){
        throw::runtime_error( "Step 1 cannot be executed: step 0 not set (starting data insertion)." );
    }

    // Generate the interleaving relative partition index arrays.
    vector< vector< unsigned int > > index_arrs = 
        LM_UTIL::gen_2_partit_idx_arr( this->fr_idx_arr.size() );
    vector< unsigned int > f1IdxVec = index_arrs.at(0);
    vector< unsigned int > f2IdxVec = index_arrs.at(1);

    // Continue the partitioning process through the standard function.
    step1_fData_partition( f1IdxVec, f2IdxVec );

}


void LM_eng::step1_fData_partition( const vector<unsigned int>& f1IdxVec, 
    const vector<unsigned int>& f2IdxVec ){
    
    if( !flag0_data_set ){
        throw::runtime_error( "Step 1 cannot be executed: step 0 not set (starting data insertion)." );
    }

    // Obtain the size of the reduced frequency set.
    unsigned int fr_len = this->fr_idx_arr.size();
    unsigned int f1_size = f1IdxVec.size();
    unsigned int f2_size = f2IdxVec.size();

    // Check if total entries between the two partitions add up to the
    // reduced f set size.
    if( f1_size + f2_size != fr_len ){
        throw std::invalid_argument( "The total number of entries from the two partitions must match the reduced f set size." );
    }

    // Create a full array from the concatenated f1 and f2 index arrays.
    vector< unsigned int > full_idx_vec = f1IdxVec;
    full_idx_vec.reserve( fr_len );
    full_idx_vec.insert( full_idx_vec.end(), f2IdxVec.begin(), f2IdxVec.end());
    // Sort array in ascending order.
    std::sort( full_idx_vec.begin(), full_idx_vec.end() );

    // Check for out of range indices.
    if( full_idx_vec.at( fr_len - 1 ) >= fr_len ){
        throw std::out_of_range( "Out of range index detected in the index arrays." );
    }
    // Check for repeated indices.
    for( unsigned int z = 1; z < fr_len - 1; z++ ){
        if( full_idx_vec[z] == full_idx_vec[z+1] ){
            throw std::invalid_argument( "Repeated indices found in the index vectors (must be non repeated in BOTH vectors)." );
        }
    }
    
    // Generate the f1 partition index array with respect to the original fdata.
    this->partit1IdxArr.clear();
    this->partit1IdxArr.reserve( f1_size );
    for ( unsigned int z = 0; z < f1IdxVec.size(); z++ ) {
        this->partit1IdxArr.push_back( fr_idx_arr[ f1IdxVec[z] ] );
    }
    // Generate the f2 partition index array with respect to the original fdata.
    this->partit2IdxArr.clear();
    this->partit2IdxArr.reserve( f2_size );
    for ( unsigned int z = 0; z < f2IdxVec.size(); z++ ) {
        this->partit2IdxArr.push_back( fr_idx_arr[ f2IdxVec[z] ] );
    }

    // Check for DC point in either partitions.
    this->f1_has_DC_pt = this->myFData.get_fval_at( this->partit1IdxArr.at(0) ) == 0;
    this->f2_has_DC_pt = this->myFData.get_fval_at( this->partit2IdxArr.at(0) ) == 0;

    // Set the tracking flag for step 1.
    this->flag1_data_prep = true;

}


void LM_eng::step2_LM_construct(){

    if( !flag1_data_prep ){
        throw::runtime_error( "Step 2 cannot be executed: step 1 not set (data prep)." );
    }

    // Create a fData partitions.
    shared_ptr<fData> myFr1 = this->myFData.red_partit( this->partit1IdxArr );
    shared_ptr<fData> myFr2 = this->myFData.red_partit( this->partit2IdxArr );

    // Generate the two partitions with their complex conjugates inserted 
    // in interleaving fashion.
    shared_ptr<fData> myFrc1 = myFr1->gen_cplx_conj_comb();
    shared_ptr<fData> myFrc2 = myFr2->gen_cplx_conj_comb();

    // Construct the Loewner Matrix using the two cconj injected partitions.
    try{
        this->LM = *LM_UTIL::build_LM( *myFrc1, *myFrc2 );
        // Construct the Loewner Matrix using the two cconj injected partitions.
        this->SLM = *LM_UTIL::build_SLM( *myFrc1, *myFrc2 );
        // Construct the W matrix vector using partition 1.
        this->W = *LM_UTIL::build_W( *myFrc1 );
        // Construct the F matrix vector using partition 2.
        this->F = *LM_UTIL::build_F( *myFrc2 );
    }catch(...){
        cerr << "step2_LM_construct exception rethrow log." << endl;
        // Rethrow exception.
        throw;
    }

    // Set the tracking flag for step 1.
    this->flag2_LM_const = true;

}


void LM_eng::step3_LM_re_trans(){
    
    if( !flag2_LM_const ){
        throw::runtime_error( "Step 3 cannot be executed: step 2 not set (LM construction)." );
    }

    // General partition data characteristics.
    unsigned int out_cnt = this->myFData.get_out_cnt();
    unsigned int in_cnt = this->myFData.get_in_cnt();
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
        throw std::runtime_error( "Step 3 has failed due to non-negligeable imaginary part remaining in transformed matrices." );
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


void LM_eng::step3skip2_LM_re_construct(){

    if( !flag1_data_prep ){
        throw::runtime_error( "Step 3 skipping step 2 cannot be executed: step 1 not set (data prep)." );
    }

    // General partition data characteristics.
    unsigned int out_cnt = this->myFData.get_out_cnt();
    unsigned int in_cnt = this->myFData.get_in_cnt();
    unsigned int fr1_len = this->partit1IdxArr.size();
    unsigned int fr2_len = this->partit2IdxArr.size();

    // Create a fData partitions.
    shared_ptr<fData> myFr1 = this->myFData.red_partit( this->partit1IdxArr );
    shared_ptr<fData> myFr2 = this->myFData.red_partit( this->partit2IdxArr );

    // Generate the two partitions with their complex conjugates inserted 
    // in interleaving fashion.
    shared_ptr<fData> myFrc1 = myFr1->gen_cplx_conj_comb();
    shared_ptr<fData> myFrc2 = myFr2->gen_cplx_conj_comb();

    // Obtain the number of complex conjugated f data points.
    unsigned int frc1_len = myFrc1->get_f_cnt();
    unsigned int frc2_len = myFrc2->get_f_cnt();

    // Obtain the expected height and width of the final real LMs.
    unsigned int LM_h = frc2_len*out_cnt;
    unsigned int LM_w = frc1_len*in_cnt;

    // Obtain purely real defintion of the matrices.
    this->LM_re = Eigen::MatrixXd( LM_h, LM_w );
    this->SLM_re = Eigen::MatrixXd( LM_h, LM_w );
    this->W_re = Eigen::MatrixXd( out_cnt, LM_w );
    this->F_re = Eigen::MatrixXd( LM_h, in_cnt );

    

    unsigned int lead_x = 0;
    unsigned int lead_y = 0;

    // Define square root of 2 that is going to be repeatedly reused.
    double sqrt_of_2 = std::sqrt(2);

    // Repeated temporary variables.
    complex<double> f2_i, f1_j;
    Eigen::MatrixXcd S1_j = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd S2_i = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd LMt_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd SLMt_ij = Eigen::MatrixXcd( out_cnt, in_cnt );

    // DC point affected portion computation.
    if( this->f1_has_DC_pt ){

        // Get the DC point data.
        Eigen::MatrixXcd f1_dc_data = myFr1->get_cplxData_at_f(0);
        
        // Set block column index to first column block.
        lead_y = 0;

        this->W_re.block( lead_x, lead_y, out_cnt, in_cnt ) = f1_dc_data.real();

        for( unsigned int i = 0; i < fr2_len; i++ ){

            // Current partition #2 frequency and data.
            f2_i = myFr2->get_cplx_f_at(i);
            S2_i = myFr2->get_cplxData_at_f(i);

            LMt_ij = S2_i - f1_dc_data;
            LMt_ij = sqrt_of_2 * LMt_ij/( f2_i );

            SLMt_ij = sqrt_of_2 * S2_i;

            lead_x = 2*( i*out_cnt );
            this->LM_re.block( lead_x, lead_y, out_cnt, in_cnt ) = LMt_ij.real();
            this->LM_re.block( lead_x + out_cnt, lead_y, out_cnt, in_cnt ) = -1*LMt_ij.imag();

            this->SLM_re.block( lead_x, lead_y, out_cnt, in_cnt ) = SLMt_ij.real();
            this->SLM_re.block( lead_x + out_cnt, lead_y, out_cnt, in_cnt ) = -1*SLMt_ij.imag();

        }


    }else if( this->f2_has_DC_pt ){

        // Get the DC point data.
        Eigen::MatrixXcd f2_dc_data = myFr2->get_cplxData_at_f(0);
        // Set block row index to first row block.
        lead_x = 0;

        this->F_re.block( lead_x, lead_y, out_cnt, in_cnt ) = f2_dc_data.real();

        for( unsigned int j = 0; j < fr1_len; j++ ){

            // Current partition #1 frequency and data.
            f1_j = myFr1->get_cplx_f_at(j);
            S1_j = myFr1->get_cplxData_at_f(j);

            LMt_ij = f2_dc_data - S1_j;
            LMt_ij = -1 * sqrt_of_2 * LMt_ij/( f1_j );

            SLMt_ij = sqrt_of_2 * S1_j;

            lead_y = 2*( j*in_cnt );
            this->LM_re.block( lead_x, lead_y, out_cnt, in_cnt ) = LMt_ij.real();
            this->LM_re.block( lead_x, lead_y + in_cnt, out_cnt, in_cnt ) = LMt_ij.imag();

            this->SLM_re.block( lead_x, lead_y, out_cnt, in_cnt ) = SLMt_ij.real();
            this->SLM_re.block( lead_x, lead_y + in_cnt, out_cnt, in_cnt ) = SLMt_ij.imag();

        }

    }


    // Remaining standard LM computations.
    
    // Define indexing offsets to take into account of the DC point.
    unsigned x_lead_offset = 0;   unsigned y_lead_offset = 0;
    unsigned int i_offset = 0;    unsigned int j_offset = 0;
    if( this->f1_has_DC_pt ){
        y_lead_offset = in_cnt;
        j_offset = 1;
    }else if( this->f2_has_DC_pt ){
        x_lead_offset = out_cnt;
        i_offset = 1;
    }

    Eigen::MatrixXcd LM_a_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd LM_b_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd SLM_a_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd SLM_b_ij = Eigen::MatrixXcd( out_cnt, in_cnt );

    for( unsigned int i = i_offset; i < fr2_len; i++ ){
        
        lead_x = 2*( ( i - i_offset )*out_cnt ) + x_lead_offset;

        for( unsigned int j = j_offset; j < fr1_len; j++ ){

            lead_y = 2*( ( j - j_offset )*in_cnt ) + y_lead_offset;

            // Current partition #2 frequency and data.
            f2_i = myFr2->get_cplx_f_at( i );
            S2_i = myFr2->get_cplxData_at_f( i );
            // Current partition #1 frequency and data.
            f1_j = myFr1->get_cplx_f_at( j );
            S1_j = myFr1->get_cplxData_at_f( j );

            // Current block LM pieces computation.
            LM_a_ij = ( S2_i - S1_j )/( f2_i - f1_j );
            LM_b_ij = ( S2_i - S1_j.conjugate() )/( f2_i - conj( f1_j ) );
            // Current block SLM pieces computation.
            SLM_a_ij = ( f2_i*S2_i - f1_j*S1_j )/( f2_i - f1_j );
            SLM_b_ij = ( f2_i*S2_i - conj( f1_j )*S1_j.conjugate() )/
                ( f2_i - conj( f1_j ) );

            // LM current block computation.
            this->LM_re.block( lead_x, lead_y, out_cnt, in_cnt ) = 
                LM_a_ij.real() + LM_b_ij.real();
            this->LM_re.block( lead_x, lead_y + in_cnt, out_cnt, in_cnt ) = 
                LM_a_ij.imag() - LM_b_ij.imag();
            this->LM_re.block( lead_x + out_cnt, lead_y, out_cnt, in_cnt ) = 
                - LM_a_ij.imag() - LM_b_ij.imag();
            this->LM_re.block( lead_x + out_cnt, lead_y + in_cnt, out_cnt, in_cnt ) = 
                LM_a_ij.real() - LM_b_ij.real();
            // SLM current block computation.
            this->SLM_re.block( lead_x, lead_y, out_cnt, in_cnt ) = 
                SLM_a_ij.real() + SLM_b_ij.real();
            this->SLM_re.block( lead_x, lead_y + in_cnt, out_cnt, in_cnt ) = 
                SLM_a_ij.imag() - SLM_b_ij.imag();
            this->SLM_re.block( lead_x + out_cnt, lead_y, out_cnt, in_cnt ) = 
                - SLM_a_ij.imag() - SLM_b_ij.imag();
            this->SLM_re.block( lead_x + out_cnt, lead_y + in_cnt, out_cnt, in_cnt ) = 
                SLM_a_ij.real() - SLM_b_ij.real();

        }
    }

    Eigen::MatrixXcd tmp = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Lead indices reset.
    lead_y = 0;    lead_x = 0;

    for( unsigned int i = i_offset; i < fr2_len; i++ ){

        lead_x = 2*( ( i - i_offset )*out_cnt ) + x_lead_offset;

        tmp = myFr2->get_cplxData_at_f(i);

        this->F_re.block( lead_x, lead_y, out_cnt, in_cnt ) = sqrt_of_2*tmp.real();
        this->F_re.block( lead_x + out_cnt, lead_y, out_cnt, in_cnt ) = -sqrt_of_2*tmp.imag();

    }

    // Lead indices reset.
    lead_y = 0;    lead_x = 0;

    for( unsigned int j = j_offset; j < fr1_len; j++ ){

        lead_y = 2*( ( j - j_offset )*in_cnt ) + y_lead_offset;

        S1_j = myFr1->get_cplxData_at_f(j);

        this->W_re.block( lead_x, lead_y, out_cnt, in_cnt ) = sqrt_of_2*S1_j.real();
        this->W_re.block( lead_x, lead_y + in_cnt, out_cnt, in_cnt ) = sqrt_of_2*S1_j.imag();

    }

    this->flag3_re_trans = true;

}


void LM_eng::step4_LM_pencil_SVD(){

    if( !flag3_re_trans ){
        throw::runtime_error( "Step 4 cannot be executed: step 3 not set (LM real transform)." );
    }

    double old_ref_f_mag = this->ref_f_mag;
    // Obtain a reference frequency value.
    this->ref_f_mag = this->myFData.get_fval_at( this->fr_idx_arr[this->fr_idx_arr.size() - 1] );
    
    try{
        step4_LM_pencil_SVD( ref_f_mag );
    }catch( ... ){
        cerr << "step4_LM_pencil_SVD exception rethrow log." << endl;
        // Revert reference frequency.
        this->ref_f_mag = old_ref_f_mag;
        throw;
    }

}


void LM_eng::step4_LM_pencil_SVD( double f_ref ){

    if( !flag3_re_trans ){
        throw::runtime_error( "Step 4 cannot be executed: step 3 not set (LM real transform)." );
    }

    shared_ptr<Eigen::MatrixXd> LM_pen;
    // Construct the LM pencil.
    try{
        LM_pen = LM_UTIL::build_LM_pencil( this->ref_f_mag, this->LM_re, this->SLM_re );
    }catch( ... ){
        cerr << "step4_LM_pencil_SVD exception rethrow log." << endl;
        throw;
    }

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
    unsigned int out_cnt = this->myFData.get_out_cnt();

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

// void LM_eng::set_fData( const fData& inData ){

//     // Set the engine's data with the given data.
//     this->myFData = inData;

//     // Reset the flags.
//     this->flag0_data_set = true;
//     this->flag1_data_prep = false;
//     this->flag2_LM_const = false;
//     this->flag3_re_trans = false;
//     this->flag4_pen_SVD = false;

// }

fData LM_eng::get_fData() const{
    if( !this->flag0_data_set ){
        throw std::runtime_error( "Cannot return fData: step1 (data set) has not been set." );
    }
    return this->myFData;
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
    return this->myFData.red_partit( this->partit1IdxArr );
}
shared_ptr<fData> LM_eng::get_Fr2() const{
    if( !this->flag1_data_prep ){
        throw std::runtime_error( "Cannot return reduced frequency partition 2: step1 (data preparation) has not been set." );
    }
    return this->myFData.red_partit( this->partit2IdxArr );
}

shared_ptr<fData> LM_eng::get_Frc1() const{
    if( !this->flag1_data_prep ){
        throw std::runtime_error( "Cannot return reduced complex conjugate frequency partition 1: step1 (data preparation) has not been set." );
    }
    return this->myFData.red_partit( this->partit1IdxArr )->gen_cplx_conj_comb();
}
shared_ptr<fData> LM_eng::get_Frc2() const{
    if( !this->flag1_data_prep ){
        throw std::runtime_error( "Cannot return reduced complex conjugate frequency partition 2: step1 (data preparation) has not been set." );
    }
    return this->myFData.red_partit( this->partit2IdxArr )->gen_cplx_conj_comb();
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


Eigen::MatrixXd LM_eng::get_LM_re() const{
    if( !this->flag3_re_trans){
        throw std::runtime_error( "Cannot return real LM: step3 (LM real transform) has not been set." );
    }
    return this->LM_re;
}
Eigen::MatrixXd LM_eng::get_SLM_re() const{
    if( !this->flag3_re_trans ){
        throw std::runtime_error( "Cannot return real SLM: step3 (LM real transform) has not been set." );
    }
    return this->SLM_re;
}
Eigen::MatrixXd LM_eng::get_W_re() const{
    if( !this->flag3_re_trans ){
        throw std::runtime_error( "Cannot return real W: step3 (LM real transform) has not been set." );
    }
    return this->W_re;
}
Eigen::MatrixXd LM_eng::get_F_re() const{
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



vector< vector< unsigned int > > LM_UTIL::gen_2_partit_idx_arr( unsigned int origSize ){
    
    // Generate the frequency partition index arrays (Interleaving).
    vector< unsigned int > tmp1 = utils::gen_even_idx_arr( 0, origSize - 1 );
    vector< unsigned int > tmp2 = utils::gen_odd_idx_arr( 0, origSize - 1 );
    // Initialize the return vector with the two index arrays.
    vector< vector< unsigned int > > retVec = { tmp1, tmp2 };

    return retVec;

}



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


shared_ptr<Eigen::MatrixXd> LM_UTIL::build_LM_re( const fData& myFr1, const fData& myFr2 ){
    
    // Obtain the size of the two partitions.
    unsigned int f1Size = myFr1.get_f_cnt();
    unsigned int f2Size = myFr2.get_f_cnt();
    if( f1Size == 0 || f2Size == 0 ){
        throw std::invalid_argument( "Empty frequency data inputs are not allowed for constructing the Loewner Matrix." );
    }

    // General partition data characteristics.
    unsigned int out_cnt = myFr1.get_out_cnt();
    unsigned int in_cnt = myFr1.get_in_cnt();
    unsigned int fr1_len = myFr1.get_f_cnt();
    unsigned int fr2_len = myFr2.get_f_cnt();
    // Determine if the DC point is present.
    bool f1_has_DC_pt = myFr1.get_fval_at(0) == 0;
    bool f2_has_DC_pt = myFr2.get_fval_at(0) == 0;

    // Obtain the number of complex conjugated f data points.
    unsigned int frc1_len = fr1_len*2;
    if( f1_has_DC_pt ){ frc1_len--; }
    unsigned int frc2_len = fr2_len*2;
    if( f2_has_DC_pt ){ frc2_len--; }

    // Obtain the expected height and width of the final real LMs.
    unsigned int LM_h = frc2_len*out_cnt;
    unsigned int LM_w = frc1_len*in_cnt;

    shared_ptr<Eigen::MatrixXd> LM_re;
    // Obtain purely real defintion of the matrices.
    *LM_re = Eigen::MatrixXd( LM_h, LM_w );

    unsigned int lead_x = 0;
    unsigned int lead_y = 0;

    // Define square root of 2 that is going to be repeatedly reused.
    double sqrt_of_2 = std::sqrt(2);
    // Repeated temporary variables.
    complex<double> f2_i, f1_j;
    Eigen::MatrixXcd S1_j = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd S2_i = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd LMt_ij = Eigen::MatrixXcd( out_cnt, in_cnt );

    // Define block matrix lead indices.
    unsigned int lead_x = 0;
    unsigned int lead_y = 0;
    
    // DC point affected portion computation.
    if( f1_has_DC_pt ){

        // Get the DC point data.
        Eigen::MatrixXcd f1_dc_data = myFr1.get_cplxData_at_f(0);
        
        // Set block column index to first column block.
        lead_y = 0;

        for( unsigned int i = 0; i < fr2_len; i++ ){

            // Current partition #2 frequency and data.
            f2_i = myFr2.get_cplx_f_at(i);
            S2_i = myFr2.get_cplxData_at_f(i);

            LMt_ij = S2_i - f1_dc_data;
            LMt_ij = sqrt_of_2 * LMt_ij/( f2_i );

            lead_x = 2*( i*out_cnt );
            LM_re->block( lead_x, lead_y, out_cnt, in_cnt ) = LMt_ij.real();
            LM_re->block( lead_x + out_cnt, lead_y, out_cnt, in_cnt ) = -1*LMt_ij.imag();

        }


    }else if( f2_has_DC_pt ){

        // Get the DC point data.
        Eigen::MatrixXcd f2_dc_data = myFr2.get_cplxData_at_f(0);
        // Set block row index to first row block.
        lead_x = 0;

        for( unsigned int j = 0; j < fr1_len; j++ ){

            // Current partition #1 frequency and data.
            f1_j = myFr1.get_cplx_f_at(j);
            S1_j = myFr1.get_cplxData_at_f(j);

            LMt_ij = f2_dc_data - S1_j;
            LMt_ij = -1 * sqrt_of_2 * LMt_ij/( f1_j );

            lead_y = 2*( j*in_cnt );
            LM_re->block( lead_x, lead_y, out_cnt, in_cnt ) = LMt_ij.real();
            LM_re->block( lead_x, lead_y + in_cnt, out_cnt, in_cnt ) = LMt_ij.imag();

        }

    }

    // Remaining standard LM computations.
    
    // Define indexing offsets to take into account of the DC point.
    unsigned x_lead_offset = 0;   unsigned y_lead_offset = 0;
    unsigned int i_offset = 0;    unsigned int j_offset = 0;
    if( f1_has_DC_pt ){
        y_lead_offset = in_cnt;
        j_offset = 1;
    }else if( f2_has_DC_pt ){
        x_lead_offset = out_cnt;
        i_offset = 1;
    }

    Eigen::MatrixXcd LM_a_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd LM_b_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    
    for( unsigned int i = i_offset; i < fr2_len; i++ ){
        
        lead_x = 2*( ( i - i_offset )*out_cnt ) + x_lead_offset;

        for( unsigned int j = j_offset; j < fr1_len; j++ ){

            lead_y = 2*( ( j - j_offset )*in_cnt ) + y_lead_offset;

            // Current partition #2 frequency and data.
            f2_i = myFr2.get_cplx_f_at( i );
            S2_i = myFr2.get_cplxData_at_f( i );
            // Current partition #1 frequency and data.
            f1_j = myFr1.get_cplx_f_at( j );
            S1_j = myFr1.get_cplxData_at_f( j );

            // Current block LM pieces computation.
            LM_a_ij = ( S2_i - S1_j )/( f2_i - f1_j );
            LM_b_ij = ( S2_i - S1_j.conjugate() )/( f2_i - conj( f1_j ) );

            // LM current block computation.
            LM_re->block( lead_x, lead_y, out_cnt, in_cnt ) = 
                LM_a_ij.real() + LM_b_ij.real();
            LM_re->block( lead_x, lead_y + in_cnt, out_cnt, in_cnt ) = 
                LM_a_ij.imag() - LM_b_ij.imag();
            LM_re->block( lead_x + out_cnt, lead_y, out_cnt, in_cnt ) = 
                - LM_a_ij.imag() - LM_b_ij.imag();
            LM_re->block( lead_x + out_cnt, lead_y + in_cnt, out_cnt, in_cnt ) = 
                LM_a_ij.real() - LM_b_ij.real();

        }
    }

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

