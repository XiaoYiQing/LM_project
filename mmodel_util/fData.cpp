#include "fData.h"



using namespace std;



// ====================================================================== >>>>>
//      Class Enum "METRIC_PREFIX" Help Functions
// ====================================================================== >>>>>

string fData::get_METRIC_PREFIX_Str( fData::METRIC_PREFIX tar_METRIC_PREFIX ){
    return string( magic_enum::enum_name( tar_METRIC_PREFIX ) );
}

fData::METRIC_PREFIX fData::get_METRIC_PREFIX_AtIdx( int idx ){
    if( idx >= 0 && idx < fData::METRIC_PREFIX_Count ){
        return static_cast<fData::METRIC_PREFIX>(idx);
    }else{
        throw std::out_of_range( "Specific enum index is out of range." );
    }
}

fData::METRIC_PREFIX fData::get_METRIC_PREFIX( string strSymbol ){
    auto tmp_enum = magic_enum::enum_cast<METRIC_PREFIX>(strSymbol);
    if( tmp_enum.has_value() ){
        return tmp_enum.value();
    }else{
        return METRIC_PREFIX::NONE;
    }
}

double fData::get_METRIC_PREFIX_val( METRIC_PREFIX tar_METRIC_PREFIX ){

    // p, n, μ, m, c, d, da, h, k, M, G, T
    double retVal = 1;
    switch( tar_METRIC_PREFIX ){
        case METRIC_PREFIX::p:
            retVal = 1e-12;  break;
        case METRIC_PREFIX::n:
            retVal = 1e-9;   break;
        case METRIC_PREFIX::μ:
            retVal = 1e-6;   break;
        case METRIC_PREFIX::mu:
            retVal = 1e-6;   break;
        case METRIC_PREFIX::m:
            retVal = 1e-3;   break;
        case METRIC_PREFIX::c:
            retVal = 1e-2;   break;
        case METRIC_PREFIX::d:
            retVal = 1e-1;   break;
        case METRIC_PREFIX::da:
            retVal = 1e1;   break;
        case METRIC_PREFIX::h:
            retVal = 1e2;   break;
        case METRIC_PREFIX::k:
            retVal = 1e3;   break;
        case METRIC_PREFIX::M:
            retVal = 1e6;   break;
        case METRIC_PREFIX::G:
            retVal = 1e9;   break;
        case METRIC_PREFIX::T:
            retVal = 1e12;   break;
        case METRIC_PREFIX::NONE:
            retVal = 1e0;   break;
        default:
            throw std::out_of_range( "Target prefix has no equivalent value defined currently by this function." );
    };

    return retVal;

}

fData::METRIC_PREFIX fData::get_METRIC_PREFIX_next( METRIC_PREFIX tar_METRIC_PREFIX, bool higher ){

    // Obtain the index equivalent to the target enum.
    auto currIdx_tmp = magic_enum::enum_index( tar_METRIC_PREFIX );
    int currIdx = -1;
    if( currIdx_tmp ){
        currIdx = *currIdx_tmp;
    }else{
        throw std::invalid_argument( "Given metric prefix is unrecognized. Abort." );
    }

    // Initialize next prefix.
    METRIC_PREFIX nextPrefix = static_cast<fData::METRIC_PREFIX>(-1);
    // Obtain the next higher metric prefix.
    try{
        if( higher ){
            nextPrefix = get_METRIC_PREFIX_AtIdx( currIdx + 1 );
        }else{
            nextPrefix = get_METRIC_PREFIX_AtIdx( currIdx - 1 );
        }
    }catch( ... ){
        throw;
    }

    return nextPrefix;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Class Enum "FDATA_TYPE" Help Functions
// ====================================================================== >>>>>

string fData::get_FDATA_TYPE_Str( fData::FDATA_TYPE tar_FDATA_TYPE ){
    return string( magic_enum::enum_name( tar_FDATA_TYPE ) );
}

fData::FDATA_TYPE fData::get_FDATA_TYPE_AtIdx( int idx ){
    if( idx >= 0 && idx < fData::FDATA_TYPE_Count ){
        return static_cast<fData::FDATA_TYPE>(idx);
    }else{
        cout << "Invalid int index for accessing enum \"FDATA_TYPE\"." << endl;
        return static_cast<fData::FDATA_TYPE>(-1);
    }
}

fData::FDATA_TYPE fData::get_FDATA_TYPE( string strSymbol ){
    auto tmp_enum = magic_enum::enum_cast<FDATA_TYPE>(strSymbol);
    if( tmp_enum.has_value() ){
        return tmp_enum.value();
    }else{
        return FDATA_TYPE::NONE;
    }
}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Class Enum "FDATA_FORMAT" Help Functions
// ====================================================================== >>>>>

string fData::get_FDATA_FORMAT_Str( fData::FDATA_FORMAT tar_FDATA_FORMAT ){
    return string( magic_enum::enum_name( tar_FDATA_FORMAT ) );
}

fData::FDATA_FORMAT fData::get_FDATA_FORMAT_AtIdx( int idx ){
    if( idx >= 0 && idx < fData::FDATA_FORMAT_Count ){
        return static_cast<fData::FDATA_FORMAT>(idx);
    }else{
        cout << "Invalid int index for accessing enum \"FDATA_FORMAT\"." << endl;
        return static_cast<fData::FDATA_FORMAT>(-1);
    }
}

fData::FDATA_FORMAT fData::get_FDATA_FORMAT( string strSymbol ){
    auto tmp_enum = magic_enum::enum_cast<FDATA_FORMAT>(strSymbol);
    if( tmp_enum.has_value() ){
        return tmp_enum.value();
    }else{
        return FDATA_FORMAT::NONE;
    }
}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

fData::fData(){

    // Native type variables initialization.
    IOcnt[0] = 0;
    IOcnt[1] = 0;

    f_pref = METRIC_PREFIX::NONE;
    fD_type = FDATA_TYPE::NONE;
    fD_format = FDATA_FORMAT::NONE;
    systImp = 50;

}


fData::fData( Eigen::VectorXd& f_vec, Matrix3DXd& Xr_vec, Matrix3DXd& Xi_vec ){

    if( f_vec.size() == 0 ){
        throw std::invalid_argument( "This constructor does not allow an empty initialization." );
    }
    if( f_vec.size() != Xr_vec.levels() || f_vec.size() != Xi_vec.levels() ){
        throw std::invalid_argument( "The frequency vector and the two data matrix vectors must have the same number of entries." );
    }
    if( !Matrix3DXd::consist_check( Xr_vec ) || !Matrix3DXd::consist_check( Xi_vec ) ){
        throw std::invalid_argument( "Both data vectors must have consistent 2D matrices." );
    }
    if( !Matrix3DXd::same_size( Xr_vec.at(0), Xi_vec.at(0) ) ){
        throw std::invalid_argument( "Size inconsistency between the two data matrix vectors." );
    }

    this->f_vec = f_vec;
    this->Xr_vec = Xr_vec;
    this->Xi_vec = Xi_vec;

    this->IOcnt[0] = Xr_vec.cols();
    this->IOcnt[1] = Xr_vec.rows();
    this->f_pref = METRIC_PREFIX::NONE;
    this->fD_type = FDATA_TYPE::NONE;
    this->fD_format = FDATA_FORMAT::NONE;
    this->systImp = 50;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Data Editing
// ====================================================================== >>>>>

void fData::reInit( unsigned int out_cnt, unsigned int in_cnt, unsigned int f_cnt ){

    if( out_cnt == 0 || in_cnt == 0 || f_cnt == 0 ){
        throw std::domain_error( "Cannot re-initialize with empty data vectors." );
    }

    this->f_vec.resize( f_cnt );
    this->f_vec.setZero();
    this->Xr_vec.reInit( out_cnt, in_cnt, f_cnt );
    this->Xi_vec.reInit( out_cnt, in_cnt, f_cnt );
    this->IOcnt[0] = in_cnt;
    this->IOcnt[1] = out_cnt;

}

void fData::copy_data( fData& tarObj, const fData& refObj ){

    tarObj.IOcnt[0] = refObj.IOcnt[0];
    tarObj.IOcnt[1] = refObj.IOcnt[1];

    tarObj.f_vec = refObj.f_vec;
    tarObj.Xr_vec = refObj.Xr_vec;
    tarObj.Xi_vec = refObj.Xi_vec;

}

void fData::copy_settings( fData& tarObj, const fData& refObj ){

    tarObj.f_pref = refObj.f_pref;
    tarObj.fD_format = refObj.fD_format;
    tarObj.fD_type = refObj.fD_type;
    tarObj.systImp = refObj.systImp;

}

void fData::data_format_Switch( FDATA_FORMAT newFormat ){

    //NONE, DB, MA, RI

    // If specified format is the same as the current one, no change.
    if( this->fD_format == newFormat ){
        return;
    }

    switch( newFormat ){

    case( FDATA_FORMAT::DB ):

        if( this->fD_format == FDATA_FORMAT::MA ){

            this->Xr_vec.elem_log10();
            Xr_vec *= 20;
            
        }else if( this->fD_format == FDATA_FORMAT::RI ){

            // Compute the data phase.
            Matrix3DXd quot = Matrix3DXd::elem_phase_comp( Xi_vec, Xr_vec );

            // Compute the decibel magnitudes.
            Xr_vec.elem_pow(2);
            Xi_vec.elem_pow(2);
            Xr_vec = Xr_vec + Xi_vec;
            Xr_vec.elem_log10();
            Xr_vec *= 10;   // Decibel data assigned to first data matrix.

            // Assigne the data phase to second data matrix.
            Xi_vec = quot;
            
        }else{
            throw std::invalid_argument( "fData::data_format_Switch: specified data format change impossible." );
        }

        break;

    case( FDATA_FORMAT::MA ):

        if( this->fD_format == FDATA_FORMAT::DB ){

            Xr_vec *= 0.05;
            Xr_vec.elem_raise_pow( 10 );

        }else if( this->fD_format == FDATA_FORMAT::RI ){

            // Compute the data phase.
            Matrix3DXd quot = Matrix3DXd::elem_phase_comp( Xi_vec, Xr_vec );

            // Compute the linear magnitudes.
            Xr_vec.elem_pow(2);
            Xi_vec.elem_pow(2);
            Xr_vec = Xr_vec + Xi_vec;
            Xr_vec.elem_pow(0.5);

            // Assigne the data phase to second data matrix.
            Xi_vec = quot;

        }else{
            throw std::invalid_argument( "fData::data_format_Switch: specified data format change impossible." );
        }

        break;

    case( FDATA_FORMAT::RI ):

        if( this->fD_format == FDATA_FORMAT::DB ){

            this->Xr_vec *= 0.05;
            this->Xr_vec.elem_raise_pow( 10 );

            Matrix3DXd Xi_vec_cos = this->Xi_vec;
            Xi_vec_cos.elem_cos();
            this->Xi_vec.elem_sin();

            this->Xi_vec = this->Xr_vec*this->Xi_vec;
            this->Xr_vec = this->Xr_vec*Xi_vec_cos;

            int lol = 0;
            
        }else if( this->fD_format == FDATA_FORMAT::MA ){

            Matrix3DXd Xi_vec_cos = this->Xi_vec;
            Xi_vec_cos.elem_cos();
            this->Xi_vec.elem_sin();

            this->Xi_vec = this->Xr_vec*this->Xi_vec;
            this->Xr_vec = this->Xr_vec*Xi_vec_cos;

        }else{
            throw std::invalid_argument( "fData::data_format_Switch: specified data format change impossible." );
        }

        break;



    // The none format makes no change.
    case( FDATA_FORMAT::NONE ):
    default:
        return;
        break;

    };

    this->fD_format = newFormat;

}

void fData::data_prefix_switch( METRIC_PREFIX newPref ){

    double curr_pref_val = get_METRIC_PREFIX_val( this->f_pref );
    double new_pref_val = get_METRIC_PREFIX_val( newPref );
    
    // Apply the rescaling to the frequency vector.
    this->f_vec *= ( curr_pref_val/new_pref_val );
    // Update the prefix.
    this->f_pref = newPref;

}


void f_normalize(  ){

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Support Functions
// ====================================================================== >>>>>

bool fData::hasDC() const{
    return this->f_vec(0) == 0;
}

shared_ptr<fData> fData::red_partit_lin( unsigned int rSize ){
    
    // Generate a linear index vector.
    vector< unsigned int > fr_idx_vec = utils::gen_lin_idx_arr( 0, this->f_vec.size(), rSize );

    // Use the main function with the created sub-indexing vector.
    return this->red_partit( fr_idx_vec );

}

shared_ptr<fData> fData::red_partit( vector< unsigned int > fr_idx_vec ){

    // Obtain the size of the frequency data set.
    unsigned int rSize = fr_idx_vec.size();

    // Create the reduced set frequency vector, real data vector, and imaginary data vector.
    Eigen::VectorXd f1_vec = Eigen::VectorXd( rSize );
    Matrix3DXd f_Xr_vec, f_Xi_vec;
    f_Xr_vec.reInit( this->Xr_vec.rows(), this->Xr_vec.cols(), rSize );
    f_Xi_vec.reInit( this->Xi_vec.rows(), this->Xi_vec.cols(), rSize );
    for( unsigned int z = 0; z < rSize; z++ ){
        f1_vec(z) = f_vec( fr_idx_vec.at(z) );
        f_Xr_vec.set( z, this->Xr_vec.at( fr_idx_vec.at(z) ) );
        f_Xi_vec.set( z, this->Xi_vec.at( fr_idx_vec.at(z) ) );
    }

    // Create the return variable.
    shared_ptr<fData> retVec = std::make_shared<fData>( f1_vec, f_Xr_vec, f_Xi_vec ) ;
    retVec->IOcnt[0] = this->IOcnt[0];
    retVec->IOcnt[1] = this->IOcnt[1];
    retVec->f_pref = this->f_pref;
    retVec->fD_type = this->fD_type;
    retVec->fD_format = this->fD_format;
    retVec->systImp = this->systImp;

    return retVec;

}


vector< shared_ptr<fData> > fData::gen_2_partit() const{

    // Initialize return fData vector.
    vector< shared_ptr<fData> > retVec;

    // Obtain the size of the frequency data set.
    unsigned int fSize = this->f_vec.size();

    // Generate the frequency partition index arrays (Interleaving).
    vector< unsigned int > f1_idx_vec = utils::gen_even_idx_arr( 0, fSize - 1 );
    vector< unsigned int > f2_idx_vec = utils::gen_odd_idx_arr( 0, fSize - 1 );
    unsigned int f1Size = f1_idx_vec.size();
    unsigned int f2Size = f2_idx_vec.size();

    // Create the partition 1 frequency vector, real data vector, and imaginary data vector.
    Eigen::VectorXd f1_vec = Eigen::VectorXd( f1Size );
    Matrix3DXd f1_Xr_vec, f1_Xi_vec;
    f1_Xr_vec.reInit( this->Xr_vec.rows(), this->Xr_vec.cols(), f1Size );
    f1_Xi_vec.reInit( this->Xi_vec.rows(), this->Xi_vec.cols(), f1Size );
    for( unsigned int z = 0; z < f1Size; z++ ){
        f1_vec(z) = f_vec( f1_idx_vec.at(z) );
        f1_Xr_vec.set( z, this->Xr_vec.at( f1_idx_vec.at(z) ) );
        f1_Xi_vec.set( z, this->Xi_vec.at( f1_idx_vec.at(z) ) );
    }

    // Create the partition 2 frequency vector, real data vector, and imaginary data vector.
    Eigen::VectorXd f2_vec = Eigen::VectorXd( f2Size );
    Matrix3DXd f2_Xr_vec, f2_Xi_vec;
    f2_Xr_vec.reInit( this->Xr_vec.rows(), this->Xr_vec.cols(), f2Size );
    f2_Xi_vec.reInit( this->Xi_vec.rows(), this->Xi_vec.cols(), f2Size );
    for( unsigned int z = 0; z < f2Size; z++ ){
        f2_vec(z) = f_vec( f2_idx_vec.at(z) );
        f2_Xr_vec.set( z, this->Xr_vec.at( f2_idx_vec.at(z) ) );
        f2_Xi_vec.set( z, this->Xi_vec.at( f2_idx_vec.at(z) ) );
    }


    // Place the two partition data objects in the return vector.
    retVec.emplace_back( std::make_shared<fData>( f1_vec, f1_Xr_vec, f1_Xi_vec ) );
    retVec.emplace_back( std::make_shared<fData>( f2_vec, f2_Xr_vec, f2_Xi_vec ) );

    for( unsigned int z = 0; z < retVec.size(); z++ ){
        retVec[z]->IOcnt[0] = this->IOcnt[0];
        retVec[z]->IOcnt[1] = this->IOcnt[1];
        retVec[z]->f_pref = this->f_pref;
        retVec[z]->fD_type = this->fD_type;
        retVec[z]->fD_format = this->fD_format;
        retVec[z]->systImp = this->systImp;
    }
    

    return retVec;

}

vector< vector< unsigned int > > fData::gen_2_partit_idx_arr() const{

    // Obtain the size of the frequency data set.
    unsigned int fSize = this->f_vec.size();
    
    // Generate the frequency partition index arrays (Interleaving).
    vector< unsigned int > tmp1 = utils::gen_even_idx_arr( 0, fSize - 1 );
    vector< unsigned int > tmp2 = utils::gen_odd_idx_arr( 0, fSize - 1 );
    // Initialize the return vector with the two index arrays.
    vector< vector< unsigned int > > retVec = { tmp1, tmp2 };

    return retVec;

}


fData fData::gen_cplx_conj_set() const{

    // Initialize the copy using the current object.
    fData cplxConj_copy = *this;

    Eigen::VectorXd tmpFVec;
    Matrix3DXd tmpXr_vec;
    Matrix3DXd tmpXi_vec;

    if( cplxConj_copy.f_vec(0) == 0 ){
        tmpFVec = cplxConj_copy.f_vec.segment( 1, cplxConj_copy.f_vec.size() - 1 );
        tmpXr_vec = cplxConj_copy.Xr_vec.segment( 1, cplxConj_copy.f_vec.size() - 1 );
        tmpXi_vec = cplxConj_copy.Xi_vec.segment( 1, cplxConj_copy.f_vec.size() - 1 );
    }else{
        tmpFVec = cplxConj_copy.f_vec;
        tmpXr_vec = cplxConj_copy.Xr_vec;
        tmpXi_vec = cplxConj_copy.Xi_vec;
    }

    switch( this->fD_format ){

    case fData::FDATA_FORMAT::DB:
    case fData::FDATA_FORMAT::MA:
    case fData::FDATA_FORMAT::RI:
        // Invert the imaginary part.
        tmpFVec *= -1;
        tmpXi_vec *= -1;
        break;

    }

    cplxConj_copy.f_vec = tmpFVec;
    cplxConj_copy.Xr_vec = tmpXr_vec;
    cplxConj_copy.Xi_vec = tmpXi_vec;

    return cplxConj_copy;

}

shared_ptr<fData> fData::gen_cplx_conj_comb() const{

    bool DC_present = this->hasDC();
    unsigned int f_cnt = this->get_f_cnt();
    unsigned int in_cnt = this->IOcnt[0];
    unsigned int out_cnt = this->IOcnt[1];
    
    // Compute the new amount of frequency entries with complex conjugates included.
    unsigned f_cnt_new = 2*f_cnt;
    if( DC_present ){
        f_cnt_new--;
    }
    
    // Initialize the frequency vector.
    Eigen::VectorXd f_vec_new( f_cnt_new );
    // Initialize the frequency data matrix vectors.
    Matrix3DXd Xr_vec_new = Matrix3DXd( out_cnt, in_cnt, f_cnt_new );
    Matrix3DXd Xi_vec_new = Matrix3DXd( out_cnt, in_cnt, f_cnt_new );
    
    // Initialize frequency index variable.
    unsigned int f_idx = 0, f_new_idx = 0;
    // Deal with the DC point insertion.
    if( DC_present ){
        f_vec_new( f_new_idx ) = this->get_fval_at( f_idx );
        Xr_vec_new.set( f_new_idx, this->get_reData_at_f( f_idx ) );
        Xi_vec_new.set( f_new_idx, this->get_imData_at_f( f_idx ) );

        f_idx++;
        f_new_idx++;
    }


    while( f_idx < f_cnt ){

        f_vec_new( f_new_idx ) = f_vec( f_idx );
        f_vec_new( f_new_idx + 1 ) = - f_vec( f_idx );

        Xr_vec_new.set( f_new_idx, this->get_reData_at_f( f_idx ) );
        Xr_vec_new.set( f_new_idx + 1, this->get_reData_at_f( f_idx ) );

        Xi_vec_new.set( f_new_idx, this->get_imData_at_f( f_idx ) );
        Xi_vec_new.set( f_new_idx + 1, -1*( this->get_imData_at_f( f_idx ) ) );

        f_idx++;
        f_new_idx += 2;

    }

    // Create the return fData object and copy all settings.
    shared_ptr<fData> retFData = make_shared<fData>( f_vec_new, Xr_vec_new, Xi_vec_new );
    fData::copy_settings( *retFData, *this );
    retFData->IOcnt[0] = this->IOcnt[0];
    retFData->IOcnt[1] = this->IOcnt[1];

    return retFData;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>

void fData::set_IO_cnt( unsigned int new_in_cnt, unsigned int new_out_cnt ){
    if( new_in_cnt == 0 || new_out_cnt == 0 ){
        throw std::invalid_argument( "Number of inputs and outputs cannot be set to 0" );
    }
    this->IOcnt[0] = new_in_cnt;
    this->IOcnt[1] = new_out_cnt;
    // Delete all existing data.
    this->f_vec.resize(0);
    this->Xr_vec.clear();
    this->Xi_vec.clear();
}

void fData::set_out_cnt( unsigned int new_out_cnt ){
    if( new_out_cnt == 0 ){
        throw std::invalid_argument( "Number of outputs cannot be set to 0" );
    }
    this->IOcnt[1] = new_out_cnt;
    // Delete all existing data.
    this->f_vec.resize(0);
    this->Xr_vec.clear();
    this->Xi_vec.clear();
    
}
unsigned int fData::get_out_cnt() const
    { return this->IOcnt[1]; }

void fData::set_in_cnt( unsigned int new_in_cnt ){
    if( new_in_cnt == 0 ){
        throw std::invalid_argument( "Number of inputs cannot be set to 0" );
    }
    this->IOcnt[0] = new_in_cnt;
    // Delete all existing data.
    this->f_vec.resize(0);
    this->Xr_vec.clear();
    this->Xi_vec.clear();
    
}
unsigned int fData::get_in_cnt() const
    { return this->IOcnt[0]; }

unsigned int fData::get_f_cnt() const
    { return this->f_vec.size(); }

string fData::get_f_scale_str() const
    { return get_METRIC_PREFIX_Str( this->f_pref ); }

double fData::get_f_scale_num() const{
    return get_METRIC_PREFIX_val( this->f_pref );
}

void fData::set_fval_at( unsigned int f_idx, double f_val ){
    if( f_idx >= this->get_f_cnt() ){
        throw std::out_of_range( "Target frequency index out of range." );
    }
    this->f_vec( f_idx ) = f_val;
}
double fData::get_fval_at( unsigned int f_idx ) const{
    return this->f_vec( f_idx );
}
void fData::set_fval_block( unsigned int lead, const Eigen::VectorXd& f_blk ){

    // Empty case, which does nothing and return.
    if( f_blk.size() == 0 ){
        return;
    }

    // Check for validity of the arguments.
    if( this->get_f_cnt() <= lead + f_blk.size() ){
        throw std::out_of_range( "The current data array's size cannot accomodate the specified data block at the specified lead point." );
    }

    // Insert the data.
    for( unsigned int z = 0; z < (unsigned int) f_blk.size(); z++ ){
        this->f_vec(z) = f_blk(z);
    }

}

std::complex<double> fData::get_cplx_f_at( unsigned int f_idx ) const{
    return std::complex<double>( 0, this->f_vec( f_idx ) );
}

void fData::set_reData_at_f( unsigned int f_idx, const Eigen::MatrixXd& new_rePart ){
    if( f_idx >= this->get_f_cnt() ){
        throw std::out_of_range( "Target frequency index out of range." );
    }
    if( new_rePart.rows() != this->IOcnt[1] || new_rePart.cols() != this->IOcnt[0] ){
        throw std::invalid_argument( "Given matrix must share the expected dimensions." );
    }
    this->Xr_vec.set( f_idx, new_rePart );
}
Eigen::MatrixXd fData::get_reData_at_f( int f_idx ) const{
    return this->Xr_vec.at( f_idx );
}

void fData::set_imData_at_f( unsigned int f_idx, const Eigen::MatrixXd& new_imPart ){
    if( f_idx >= this->get_f_cnt() ){
        throw std::out_of_range( "Target frequency index out of range." );
    }
    if( new_imPart.rows() != this->IOcnt[1] || new_imPart.cols() != this->IOcnt[0] ){
        throw std::invalid_argument( "Given matrix must share the expected dimensions." );
    }
    this->Xi_vec.set( f_idx, new_imPart );
}
Eigen::MatrixXd fData::get_imData_at_f( unsigned int f_idx ) const{
    return this->Xi_vec.at( f_idx );
}


void fData::set_cplxData_at_f( unsigned int f_idx, Eigen::MatrixXcd& new_mat ){
    if( f_idx >= this->get_f_cnt() ){
        throw std::out_of_range( "Target frequency index out of range." );
    }
    if( new_mat.rows() != this->IOcnt[1] || new_mat.cols() != this->IOcnt[0] ){
        throw std::invalid_argument( "Given matrix must share the expected dimensions." );
    }
    this->Xr_vec.set( f_idx, new_mat.real() );
    this->Xi_vec.set( f_idx, new_mat.imag() );
}
Eigen::MatrixXcd fData::get_cplxData_at_f( unsigned int f_idx ) const{

    Eigen::MatrixXcd cplxMat( Xr_vec.rows(), Xr_vec.cols() );
    for( unsigned int i = 0; i < Xr_vec.rows(); i++ ) {
        for( unsigned int j = 0; j < Xr_vec.cols(); j++ ) {
            cplxMat(i, j) = std::complex<double>( Xr_vec.at(f_idx)(i, j), Xi_vec.at(f_idx)(i, j) );
        }
    }
    return cplxMat;

}

void fData::set_cplxData_block( unsigned int lead, const Matrix3DXd& newBlk_re, const Matrix3DXd& newBlk_im ){

    // Empty case, which does nothing and return.
    if( newBlk_re.isEmpty() || newBlk_im.isEmpty() ){
        return;
    }

    // Check for validity of the arguments.
    if( !Matrix3DXd::same_size( newBlk_re.at(0), newBlk_im.at(0) ) ){
        throw std::invalid_argument( "The real and imaginary part must have consistent matrix sizes." );
    }
    if( this->get_f_cnt() <= lead + newBlk_re.levels() ){
        throw std::out_of_range( "The current data array's size cannot accomodate the specified data block at the specified lead point." );
    }

    // Insert the data.
    for( unsigned int z = 0; z < newBlk_re.levels(); z++ ){
        this->Xr_vec.set( z, newBlk_re.at(z) );
        this->Xi_vec.set( z, newBlk_im.at(z) );
    }

}

Eigen::VectorXd fData::getF_vec() const{
    // Eigen::VectorXd tmp = this->f_vec;
    return this->f_vec;
}
Matrix3DXd fData::getXr_vec() const
    { return this->Xr_vec; }
Matrix3DXd fData::getXi_vec() const
    { return this->Xi_vec; }

// ====================================================================== <<<<<



// ====================================================================== >>>>>
//      Printing
// ====================================================================== >>>>>

void fData::print_to( const string& fileDir, const string& fileStem, int options ){

// ---------------------------------------------------------------------- >>>>>
//      File Name Editing
// ---------------------------------------------------------------------- >>>>>

    string fileExt = ".txt";
    string fileName = fileStem + fileExt;
    string fullFileName = fileDir + "/" + fileName;

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Stream Prep
// ---------------------------------------------------------------------- >>>>>

    // Open the file stream.
    std::ofstream file(fullFileName);
    if (!file.is_open()) {
        throw::invalid_argument( "Cannot open stream for specified file. ABORT." );
    }

    // Obtain the data count of the current frequency data.
    unsigned int data_cnt = this->f_vec.size();
    // Obtain the number of inputs and outputs.
    unsigned int row_cnt = this->get_out_cnt();
    unsigned int col_cnt = this->get_in_cnt();

    // Initialize temporary complex variable.
    complex<double> tmp_cplx = 0;

    // Set the precision of the number being printed.
    unsigned int precision = 10;
    
// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Write Info Lines
// ---------------------------------------------------------------------- >>>>>

    string fType = get_FDATA_TYPE_Str( this->fD_type );
    file << "Freq" << " ";
    for (unsigned int i = 0; i < row_cnt; ++i) {
    for (unsigned int j = 0; j < col_cnt; ++j) {
        file << fType << std::to_string(i+1) << std::to_string(j+1) << "a ";
        file << fType << std::to_string(i+1) << std::to_string(j+1) << "b ";
    }
    }
    string a_str = "NONE";
    string b_str = "NONE";
    switch( this->fD_format ){
    case FDATA_FORMAT::RI:
        a_str = "re";
        b_str = "im";
        break;
    case FDATA_FORMAT::MA:
        a_str = "mag";
        b_str = "phase(rad)";
        break;
    case FDATA_FORMAT::DB:
        a_str = "mag(dB)";
        b_str = "phase(rad)";
        break;
    }
    file << "a=[" << a_str <<"],b=[" << b_str <<"]";
    file << "\n";


// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Write Data Lines
// ---------------------------------------------------------------------- >>>>>

    for( unsigned int z = 0; z < data_cnt; z++ ){

        Eigen::MatrixXd Xr_z = this->Xr_vec.at(z);
        Eigen::MatrixXd Xi_z = this->Xi_vec.at(z);

        // switch( this->fD_format ){
        //     case FDATA_FORMAT::RI:
        //     case FDATA_FORMAT::MA:
        //     case FDATA_FORMAT::DB:
        // }

        // Print the frequency value of the current point.
        file << std::scientific << std::setprecision(precision) << this->f_vec(z) << " ";

        // Write matrix data
        for (unsigned int i = 0; i < row_cnt; ++i) {
            for (unsigned int j = 0; j < col_cnt; ++j) {
                if( Xr_z(i,j) >= 0 ){
                    file << "+";
                }
                file << std::scientific << std::setprecision(precision) << Xr_z(i,j) << " ";
                if( Xi_z(i,j) >= 0 ){
                    file << "+";
                }
                file << std::scientific << std::setprecision(precision) << Xi_z(i,j);
                if (j < col_cnt - 1) file << " ";
            }
            if (i < row_cnt - 1) file << " ";
        }
        file << "\n";

    }

// ---------------------------------------------------------------------- <<<<<

    // Close the file once all is done.
    file.close();

}

// ====================================================================== <<<<<



// ====================================================================== >>>>>
//      Data From File Utility
// ====================================================================== >>>>>

void fData::read_sXp_file( fData& tarFData, const string& fullFileName ){

// ---------------------------------------------------------------------- >>>>>
//      Full Name Parsing
// ---------------------------------------------------------------------- >>>>>

    std::filesystem::path fullFilePath(fullFileName);
    std::string fileDir = fullFilePath.parent_path().string();
    std::string fileStem = fullFilePath.stem().string();
    std::string fileExt = fullFilePath.extension().string();
    
// ---------------------------------------------------------------------- <<<<<

// ---------------------------------------------------------------------- >>>>>
//      Port Count Regex Determination
// ---------------------------------------------------------------------- >>>>>

    /*
    Regex for determining the positive integer value X in the pattern ".sXp".
    */
    regex pattern(R"(\.s(\d+)p)");
    // The match result variable.
    smatch matches;
    // The exact number of the match in string.
    string xValue;
    
    if ( regex_match( fileExt, matches, pattern ) ) {
        if ( matches.size() > 1 ) { // Check if we have a match for the integer
            xValue = matches[1]; // Get the captured group (X)
            // cout << "X in .sXp is: " << xValue << endl; // Output the value of X
        }
    } else {
        throw std::invalid_argument( "File extension is not that of a touchstone file. Parsing aborted." );
    }
    
    
    int port_cnt_tmp = -1;
    try {
        port_cnt_tmp = std::stoi( xValue );
    } catch (const std::invalid_argument& e) {
        cout << "Invalid input: The string does not contain a valid integer." << std::endl;
        cout << e.what() << endl;
        return;
    } catch (const std::out_of_range& e) {
        cout << "Invalid input: The integer is out of range." << std::endl;
        cout << e.what() << endl;
        return;
    }
    
    /*
    .s2p file case.
    The 2-port case is parsed differently from all other cases.
    */
    if( port_cnt_tmp == 2 ){
        read_s2p_file( tarFData, fullFileName );
        return;
    }

    // UPdate the port count.
    tarFData.IOcnt[0] = port_cnt_tmp;
    tarFData.IOcnt[1] = port_cnt_tmp;

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      File Parameters Read
// ---------------------------------------------------------------------- >>>>>

    // Open the input file stream.
    std::ifstream inputFile( fullFilePath );
    if( !inputFile ){
        throw std::invalid_argument( "Failed to open stream for .sXp file. Abort parsing." );
    }else{
        cout << "File \"" + fileStem + fileExt + "\": stream opened successful." << endl;
    }


    // If the first character of a line is == comm_mark, the line is a comment.
    string comm_mark = "!";
    // If the first character of a line is == opt_line_mark, the line holds the various
    // options of the data format.
    string opt_line_mark = "#";
    // The variables holding the line and word currently read, respectively.
    string line, word;

    // Define the data options' default values.
    vector<string> options = { "GHZ", "S", "MA", "R", "50" };
    
    // Boolean flag for indicating whether the line stream has reached the data lines.
    bool data_reached = false;

    while( !data_reached && getline( inputFile, line ) ){

        // Set the stream for the current line.
        istringstream iss(line);
        // Read the first word.
        iss >> word;

        // Check for comment line mark.
        if( word == comm_mark ){
            //cout << "This is a comment!" << endl;

        // Check for option line mark.
        }else if( word == opt_line_mark ){
            //cout << "This is an option!" << endl;
            int opt_idx = 0;
            while( iss >> word ){
                options.at(opt_idx) = word;
                opt_idx++;
            }

        // Check for end of consecutive series of comment and option lines.
        }else{

            data_reached = true;

        }

    }

    // If we reached the end of the file without reaching any data line, abort.
    if( !data_reached ){
        throw std::runtime_error( "The entire file has been read without reaching a data line." );
    }


    // Obtain the frequency metric prefix.
    fData::METRIC_PREFIX f_pref = METRIC_PREFIX::NONE;
    if( options.at(0).size() == 3 ){
        char keyChar = options.at(0)[0];
        fData::METRIC_PREFIX f_pref =  fData::get_METRIC_PREFIX( string( 1, keyChar ) );
        tarFData.f_pref = f_pref;
    }
    
    // Obtain the f data type.
    tarFData.fD_type = fData::get_FDATA_TYPE( options.at(1) );
    // Obtain the f data format.
    tarFData.fD_format = fData::get_FDATA_FORMAT( options.at(2) );
    // Obtain the input impedance of the measurement.
    if( options.at(3) == "R" ){
        tarFData.systImp = std::stod( options.at(4) );
    }

    // Add multiplier modifiers for the data depending on the format of the data
    // being parsed.
    double myPI = 2*std::asin(1.0);
    double Xi_fac = 1.0;
    if( tarFData.fD_format == fData::FDATA_FORMAT::DB || 
        tarFData.fD_format == fData::FDATA_FORMAT::MA ){
        Xi_fac = myPI/180;
    }

// ---------------------------------------------------------------------- <<<<<

    

// ---------------------------------------------------------------------- >>>>>
//      File Data Read
// ---------------------------------------------------------------------- >>>>>

    // Total number of parameter within the data matrix 
    // (For example, 2 by 2 S-parameters matrix has 4 individual S-parameters).
    unsigned int mat_ent_cnt = tarFData.IOcnt[0]*tarFData.IOcnt[1];

    // File parsing control variables.
    unsigned int f_idx = 0;
    unsigned int data_idx = 0;
    unsigned int res_blk_size = 200;
    unsigned int curr_vec_size = 0;
    // Temporary data value to be used during translation from string to double.
    double tmp_val = 0;

    // The frequency vector.
    vector< double > f_vec;
    // Initialize the vector of matrices.
    tarFData.Xr_vec.reInit( tarFData.IOcnt[1], tarFData.IOcnt[0], res_blk_size );
    tarFData.Xi_vec.reInit( tarFData.IOcnt[1], tarFData.IOcnt[0], res_blk_size );

    f_vec.reserve( res_blk_size );

    // Update vector size.
    curr_vec_size = tarFData.Xr_vec.levels();

    do{

        // Set the stream for the current line.
        istringstream iss( line );

        // Read the first word (frequency or comment mark).
        iss >> word;

        // Line skip if comment line.
        if( word == comm_mark ){
            continue;
        }

        // Translate the word into a double value freq.
        tmp_val = std::stod( word );
        // Save the frequency value.
        f_vec.push_back( tmp_val );

        // Read the remaining data from the first line at the current f.
        for( int j = 0; j < tarFData.IOcnt[0]; j++ ){
            // Read the next data mag.
            iss >> word;    tmp_val = std::stod( word );
            tarFData.Xr_vec.set( 0, j, f_idx, tmp_val );
            // Read the next data phase.
            iss >> word;    tmp_val = std::stod( word );
            tarFData.Xi_vec.set( 0, j, f_idx, tmp_val*Xi_fac );
        }

        // Read the remaining data lines at the current f.
        for( int i = 1; i < tarFData.IOcnt[1]; i++ ){

            // Set the stream for the next data line at the current freq. point.
            getline( inputFile, line );
            istringstream iss2( line );

            for( int j = 0; j < tarFData.IOcnt[0]; j++ ){

                // Read the next data mag.
                iss2 >> word;    tmp_val = std::stod( word );
                tarFData.Xr_vec.set( i, j, f_idx, tmp_val );
                // Read the next data phase.
                iss2 >> word;    tmp_val = std::stod( word );
                tarFData.Xi_vec.set( i, j, f_idx, tmp_val*Xi_fac );

            }
        }


        f_idx++;
        if( f_idx >= curr_vec_size ){

            f_vec.reserve( f_idx + res_blk_size );

            tarFData.Xr_vec.reserve( f_idx + res_blk_size );
            tarFData.Xi_vec.reserve( f_idx + res_blk_size );

            curr_vec_size += res_blk_size;

        }

    }while( getline( inputFile, line ) );

    // Deallocate unused reserved memory from the vectors.
    f_vec.shrink_to_fit();
    tarFData.f_vec = Eigen::Map<Eigen::VectorXd>( f_vec.data(), f_vec.size() );
    tarFData.Xr_vec.resize( f_idx );
    tarFData.Xr_vec.shrink_to_fit();
    tarFData.Xi_vec.resize( f_idx );
    tarFData.Xi_vec.shrink_to_fit();

// ---------------------------------------------------------------------- <<<<<

    cout << "File \"" + fileStem + fileExt + "\": successfully read." << endl;
    cout << "Frequency count: " << tarFData.get_f_cnt() << endl;
    cout << "Port count: " << tarFData.get_f_cnt() << endl;
    cout << endl;

    return;
    
}


void fData::read_s2p_file( fData& tarFData, const string& fullFileName ){

// ---------------------------------------------------------------------- >>>>>
//      Full Name Parsing
// ---------------------------------------------------------------------- >>>>>

    std::filesystem::path fullFilePath(fullFileName);
    std::string fileDir = fullFilePath.parent_path().string();
    std::string fileStem = fullFilePath.stem().string();
    std::string fileExt = fullFilePath.extension().string();
    
// ---------------------------------------------------------------------- <<<<<

    // UPdate the port count.
    tarFData.IOcnt[0] = 2;
    tarFData.IOcnt[1] = 2;

// ---------------------------------------------------------------------- >>>>>
//      File Parameters Read
// ---------------------------------------------------------------------- >>>>>

    // Open the input file stream.
    std::ifstream inputFile( fullFilePath );
    if( !inputFile ){
        cout << "Failed to read file." << endl;
        return;
    }else{
        cout << "File \"" + fileStem + fileExt + "\": stream opened successful." << endl;
    }


    // If the first character of a line is == comm_mark, the line is a comment.
    string comm_mark = "!";
    // If the first character of a line is == opt_line_mark, the line holds the various
    // options of the data format.
    string opt_line_mark = "#";
    // The variables holding the line and word currently read, respectively.
    string line, word;

    // Define the data options' default values.
    vector<string> options = { "GHZ", "S", "MA", "R", "50" };
    
    // Boolean flag for indicating whether the line stream has reached the data lines.
    bool data_reached = false;

    while( !data_reached && getline( inputFile, line ) ){

        // Set the stream for the current line.
        istringstream iss(line);
        // Read the first word.
        iss >> word;

        // Check for comment line mark.
        if( word == comm_mark ){
            //cout << "This is a comment!" << endl;

        // Check for option line mark.
        }else if( word == opt_line_mark ){
            //cout << "This is an option!" << endl;
            int opt_idx = 0;
            while( iss >> word ){
                options.at(opt_idx) = word;
                opt_idx++;
            }

        // Check for end of consecutive series of comment and option lines.
        }else{

            data_reached = true;

        }

    }

    // If we reached the end of the file without reaching any data line, abort.
    if( !data_reached ){
        cout << "The entire file has been read without reaching a data line." << endl;
        return;
    }


    // Obtain the frequency metric prefix.
    fData::METRIC_PREFIX f_pref = METRIC_PREFIX::NONE;
    if( options.at(0).size() == 3 ){
        char keyChar = options.at(0)[0];
        fData::METRIC_PREFIX f_pref =  fData::get_METRIC_PREFIX( string( 1, keyChar ) );
        cout << fData::get_METRIC_PREFIX_Str( f_pref ) << endl;
        tarFData.f_pref = f_pref;
    }
    
    // Obtain the f data type.
    tarFData.fD_type = fData::get_FDATA_TYPE( options.at(1) );
    // Obtain the f data format.
    tarFData.fD_format = fData::get_FDATA_FORMAT( options.at(2) );
    // Obtain the input impedance of the measurement.
    if( options.at(3) == "R" ){
        tarFData.systImp = std::stod( options.at(4) );
    }

    // Add multiplier modifiers for the data depending on the format of the data
    // being parsed.
    double myPI = 2*std::asin(1.0);
    double Xi_fac = 1.0;
    if( tarFData.fD_format == fData::FDATA_FORMAT::DB || 
        tarFData.fD_format == fData::FDATA_FORMAT::MA ){
        Xi_fac = myPI/180;
    }

// ---------------------------------------------------------------------- <<<<<

    

// ---------------------------------------------------------------------- >>>>>
//      File Data Read
// ---------------------------------------------------------------------- >>>>>

    // Total number of parameter within the data matrix 
    // (For example, 2 by 2 S-parameters matrix has 4 individual S-parameters).
    unsigned int mat_ent_cnt = tarFData.IOcnt[0]*tarFData.IOcnt[1];

    // File parsing control variables.
    unsigned int line_idx = 0;
    unsigned int data_idx = 0;
    unsigned int res_blk_size = 200;
    unsigned int curr_vec_size = 0;
    // Temporary data value to be used during translation from string to double.
    double tmp_val = 0;

    // The frequency vector.
    vector< double > f_vec;
    // Initialize the vector of matrices.
    tarFData.Xr_vec.reInit( tarFData.IOcnt[1], tarFData.IOcnt[0], res_blk_size );
    tarFData.Xi_vec.reInit( tarFData.IOcnt[1], tarFData.IOcnt[0], res_blk_size );

    f_vec.reserve( res_blk_size );

    // Update vector size.
    curr_vec_size = tarFData.Xr_vec.levels();

    do{

        // Set the stream for the current line.
        istringstream iss( line );

        // Read the frequency word.
        iss >> word;
        // Translate the word into a double value freq.
        tmp_val = std::stod( word );
        // Save the frequency value.
        f_vec.push_back( tmp_val );

        
        for( int i = 0; i < tarFData.IOcnt[1]; i++ ){
            for( int j = 0; j < tarFData.IOcnt[0]; j++ ){
                // Read the next data mag.
                iss >> word;    tmp_val = std::stod( word );
                tarFData.Xr_vec.set( i, j, line_idx, tmp_val );
                // Read the next data phase.
                iss >> word;    tmp_val = std::stod( word );
                tarFData.Xi_vec.set( i, j, line_idx, tmp_val*Xi_fac );
            }
        }

        line_idx++;
        if( line_idx >= curr_vec_size ){

            f_vec.reserve( line_idx + res_blk_size );

            tarFData.Xr_vec.reserve( line_idx + res_blk_size );
            tarFData.Xi_vec.reserve( line_idx + res_blk_size );

            curr_vec_size += res_blk_size;

        }

    }while( getline( inputFile, line ) );

    // Deallocate unused reserved memory from the vectors.
    f_vec.shrink_to_fit();
    tarFData.f_vec = Eigen::Map<Eigen::VectorXd>( f_vec.data(), f_vec.size() );
    tarFData.Xr_vec.resize( line_idx );
    tarFData.Xr_vec.shrink_to_fit();
    tarFData.Xi_vec.resize( line_idx );
    tarFData.Xi_vec.shrink_to_fit();

// ---------------------------------------------------------------------- <<<<<

    cout << "File \"" + fileStem + fileExt + "\": successfully read." << endl;
    cout << endl;
    
    return;

}


void fData::read_LTspice_Sp_file( fData& tarFData, const string& fullFileName ){
    
// ---------------------------------------------------------------------- >>>>>
//      Full Name Parsing
// ---------------------------------------------------------------------- >>>>>
    
    std::filesystem::path fullFilePath(fullFileName);
    std::string fileDir = fullFilePath.parent_path().string();
    std::string fileStem = fullFilePath.stem().string();
    std::string fileExt = fullFilePath.extension().string();

    if( fileExt != ".txt" ){
        throw std::invalid_argument( "Expected data file to be in .txt format." );
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      File First Line Read (Parameters Line)
// ---------------------------------------------------------------------- >>>>>

    // Open the input file stream.
    std::ifstream inputFile( fullFilePath );
    if( !inputFile ){
        throw std::invalid_argument( "LTspice S-parameter data file parsing failed: failed to open stream." );
    }else{
        cout << "File read successful." << endl;
    }
    
    // The variables holding the line and word currently read, respectively.
    string line, word;

    // Obtain the first line.
    getline( inputFile, line );

    // Set the stream for the current line.
    istringstream iss(line);

    // Define array of index pairs for indicating the order of S-parameter data
    // being presented in the data text file.
    vector< std::pair<unsigned int, unsigned int> > Sp_order;

    // Read this first word, which is expected to be "Freq."
    iss >> word;
    if( word != "Freq." ){
        throw::invalid_argument( "LTspice S-parameter data file has unexpected pattern. Parsing aborted." );
    }

    // regex pattern for identifying S-parameter indices;
    regex pat_Sname(R"(S\d\d\([^\)]*\))");
    // The match result variable.
    smatch matches;
    // The exact number of the match in string.
    string xValue;
    // Define the sub-string starting index and length.
    int start_index = 1;    int length = 2;

    // Initialize temporary S-parameter number.
    std::pair<unsigned int,unsigned int> tmp_S_val = 
        std::pair<unsigned int,unsigned int>(0,0);

    // Read the S-parameter labels following the "Freq. keyword."
    while( iss >> word ){
        
        if ( regex_match( word, matches, pat_Sname ) ) {
            if ( matches.size() == 1 ) { // Check if we have a match for the integer
                try {
                    xValue = matches[0]; // Get the captured group (X)
                    tmp_S_val.first = stoi( xValue.substr( start_index, 1 ) );
                    tmp_S_val.second = stoi( xValue.substr( start_index+1, 1 ) );
                    Sp_order.push_back( tmp_S_val );
                } catch (const std::invalid_argument& e) {
                    cout << "Invalid input: The string does not contain a valid integer." << std::endl;
                    cout << e.what() << endl;
                    return;
                } catch (const std::out_of_range& e) {
                    cout << "Invalid input: The integer is out of range." << std::endl;
                    cout << e.what() << endl;
                    return;
                }
            }
        } else {
            throw::invalid_argument( "LTspice S-parameter data file failed: unrecognized pattern amongst the S-param labels." );
        }

    }


    // Obtain the number of S-parameters.
    unsigned int Sp_cnt = Sp_order.size();
    double tmp_d_sqrt = sqrt( (double) Sp_cnt );
    int port_cnt_tmp = (int) tmp_d_sqrt;
    // Make sure the number of parameters is a square number.
    if( tmp_d_sqrt - port_cnt_tmp != 0 ){
        throw std::invalid_argument( "LTspice S-parameter data file failed: Number of total S-parameters is not a square number." );
    }

    // UPdate the port count.
    tarFData.IOcnt[0] = port_cnt_tmp;
    tarFData.IOcnt[1] = port_cnt_tmp;


    fData::METRIC_PREFIX f_pref = METRIC_PREFIX::NONE;
    tarFData.f_pref = f_pref;
    // Set freq data type to S by default.
    tarFData.fD_type = fData::FDATA_TYPE::S;
    // Set freq data format to dB by default.
    tarFData.fD_format = fData::FDATA_FORMAT::DB;
    // Set input impedance to 50 ohms by default.
    // NOTE: LTspice data file does not provide input impedance.
    tarFData.systImp = 50.0;
    
    // Add multiplier modifiers for the data depending on the format of the data
    // being parsed.
    double myPI = 2*std::asin(1.0);
    double Xi_fac = 1.0;
    if( tarFData.fD_format == fData::FDATA_FORMAT::DB || 
        tarFData.fD_format == fData::FDATA_FORMAT::MA ){
        Xi_fac = myPI/180;
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      File Data Read
// ---------------------------------------------------------------------- >>>>>
    
    // File parsing control variables.
    unsigned int line_idx = 0;
    unsigned int res_blk_size = 200;
    unsigned int curr_vec_size = 0;

    // The frequency vector.
    vector< double > f_vec;
    f_vec.reserve( res_blk_size );
    // Initialize the vector of matrices.
    tarFData.Xr_vec.reInit( tarFData.IOcnt[1], tarFData.IOcnt[0], res_blk_size );
    tarFData.Xi_vec.reInit( tarFData.IOcnt[1], tarFData.IOcnt[0], res_blk_size );

    // Update vector size.
    curr_vec_size = tarFData.Xr_vec.levels();

    

    // Temporary variable for storing the current double being read.
    double tmp_d = 0.0;
    string tmp_s = "";
    
    // Variable for expected magnitude and phase units (dB,degree).
    pair<string,string> exp_units = pair<string,string>( "dB", "\u00B0" );

    // Regex for filtering the magnitude value.
    regex pat_mag(R"(.*\((.*),.*)");
    // Regex for separating the magnitude numeric value from its unit.
    regex pat_mag_unit( R"(^([^a-zA-Z]*e[^a-zA-Z]*)(.*))" );
    // Regex for filtering the phase value.
    regex pat_phase(R"(.*,(.*)\))");
    // Regex for separating the phase numeric value from its unit.
    regex pat_phase_unit( R"(^([0-9.+-]+e[0-9.+-]+)(.*))" );

    
    // Reach the second line, which is the first data line.
    getline( inputFile, line );
    do{

        // Set the stream for the current line.
        istringstream iss(line);

        // Read the current frequency word.
        iss >> word;
        // Translate the word into a double value freq.
        tmp_d = std::stod( word );
        f_vec.push_back( tmp_d );

        for( unsigned int z = 0; z < Sp_cnt; z++ ){

            // Obtain the next S-parameter value in string.
            iss >> word;

            // Parse the magnitude.
            if ( regex_match( word, matches, pat_mag ) ){

                // Obtain the full magnitude string.
                xValue = matches[1];
                // Separate the magnitude numerical value from its unit.
                regex_match( xValue, matches, pat_mag_unit );
                if( matches.size() != 3 ){
                    throw std::invalid_argument( "LTspice parser encountered unexpected data pattern. ABORT." );
                }
                // Save the magnitude.
                tarFData.Xr_vec.set( Sp_order.at(z).first - 1, Sp_order.at(z).second - 1, 
                    line_idx, std::stod( matches[1] ) );
                // Check unit consistency.
                if( exp_units.first != matches[2] ){
                    throw std::invalid_argument( "Unexpected magnitude unit. Should be " + exp_units.first );
                }

            }

            // Parse the phase.
            if ( regex_match( word, matches, pat_phase ) ){

                // Obtain the full phase string.
                xValue = matches[1];
                // Separate the phase numerical value from its unit.
                regex_match( xValue, matches, pat_phase_unit );
                if( matches.size() != 3 ){
                    throw std::invalid_argument( "LTspice parser encountered unexpected data pattern. ABORT." );
                }
                // Save the phase.
                tarFData.Xi_vec.set( Sp_order.at(z).first - 1, Sp_order.at(z).second - 1, 
                    line_idx, Xi_fac*std::stod( matches[1] ) );
                // Check unit consistency.
                if( exp_units.second != matches[2] ){
                    throw std::invalid_argument( "Unexpected phase unit. Should be " + exp_units.first );
                }

            }

        }

        line_idx++;
        // Increase the data arrays reserved space if current limit reached.
        if( line_idx >= curr_vec_size ){

            f_vec.reserve( line_idx + res_blk_size );
            tarFData.Xr_vec.reserve( line_idx + res_blk_size );
            tarFData.Xi_vec.reserve( line_idx + res_blk_size );

            curr_vec_size += res_blk_size;

        }

    }while( getline( inputFile, line ) );

    // Deallocate unused reserved memory from the vectors.
    f_vec.shrink_to_fit();
    tarFData.f_vec = Eigen::Map<Eigen::VectorXd>( f_vec.data(), f_vec.size() );
    tarFData.Xr_vec.resize( line_idx );
    tarFData.Xr_vec.shrink_to_fit();
    tarFData.Xi_vec.resize( line_idx );
    tarFData.Xi_vec.shrink_to_fit();

// ---------------------------------------------------------------------- <<<<<

    return;

}

// ====================================================================== <<<<<