#include "fData.h"



using namespace std;



// ====================================================================== >>>>>
//      Class Enum "CHK_PIECE" Help Functions
// ====================================================================== >>>>>

string fData::get_METRIC_PREFIX_Str( fData::METRIC_PREFIX tar_METRIC_PREFIX ){
    return string( magic_enum::enum_name( tar_METRIC_PREFIX ) );
}

fData::METRIC_PREFIX fData::get_METRIC_PREFIX_AtIdx( int idx ){
    if( idx >= 0 && idx < fData::METRIC_PREFIX_Count ){
        return static_cast<fData::METRIC_PREFIX>(idx);
    }else{
        cout << "Invalid int index for accessing enum \"METRIC_PREFIX\"." << endl;
        return static_cast<fData::METRIC_PREFIX>(-1);
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
            retVal = 1e0;   break;
    };

    return retVal;

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

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Data Editing
// ====================================================================== >>>>>

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
            cerr << "An impossible outcome has been reached. Abort" << endl;
            return;
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
            cerr << "An impossible outcome has been reached. Abort" << endl;
            return;
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
            cerr << "An impossible outcome has been reached. Abort" << endl;
            return;
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

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Support Functions
// ====================================================================== >>>>>

    fData fData::gen_cplx_conj_set(){

        // Initialize the copy using the current object.
        fData cplxConj_copy = *this;

        

        cout << cplxConj_copy.f_vec(0) << ", " << cplxConj_copy.f_vec(1) << endl;

        Eigen::VectorXd tmpVec;
        if( cplxConj_copy.f_vec(0) == 0 ){
            tmpVec = cplxConj_copy.f_vec.segment( 1, cplxConj_copy.f_vec.size() - 1 );
        }else{
            tmpVec = cplxConj_copy.f_vec.segment( 1, cplxConj_copy.f_vec.size() - 1 );
        }

        cout << tmpVec(0) << ", " << tmpVec(1) << endl;

        cplxConj_copy.f_vec = tmpVec;

        cout << cplxConj_copy.f_vec(0) << ", " << cplxConj_copy.f_vec(1) << endl;

        int lol = 0;

        return cplxConj_copy;

    }

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>

int fData::get_out_cnt() const
    { return this->IOcnt[0]; }
int fData::get_in_cnt() const
    { return this->IOcnt[1]; }

int fData::get_f_cnt() const
    { return this->f_vec.size(); }

string fData::get_f_scale_str() const
    { return get_METRIC_PREFIX_Str( this->f_pref ); }

double fData::get_f_scale_num() const{
    return get_METRIC_PREFIX_val( this->f_pref );
}

Eigen::MatrixXd fData::get_reData_at_f( int f_idx ) const{
    return this->Xr_vec.at( f_idx );
}
Eigen::MatrixXd fData::get_imData_at_f( int f_idx ) const{
    return this->Xi_vec.at( f_idx );
}

// ====================================================================== <<<<<



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
                cout << "X: " << xValue << endl; // Output the value of X
            }
        } else {
            cout << "No match found." << endl;
            return;
        }
        
        
        int port_cnt_tmp = -1;
        try {
            port_cnt_tmp = std::stoi( xValue );
            const int loool = std::stoi( xValue );
            cout << "The integer is: " << port_cnt_tmp << std::endl;
        } catch (const std::invalid_argument& e) {
            cout << "Invalid input: The string does not contain a valid integer." << std::endl;
            cout << e.what() << endl;
            return;
        } catch (const std::out_of_range& e) {
            cout << "Invalid input: The integer is out of range." << std::endl;
            cout << e.what() << endl;
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
            cout << "Failed to read file." << endl;
            return;
        }else{
            cout << "File read successful." << endl;
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
                cout << "This is a comment!" << endl;
    
            // Check for option line mark.
            }else if( word == opt_line_mark ){
                cout << "This is an option!" << endl;
                int opt_idx = 0;
                while( iss >> word ){
                    options.at(opt_idx) = word;
                    opt_idx++;
                }
                cout << endl;
    
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
        }
        tarFData.f_pref = f_pref;
        // Obtain the f data type.
        tarFData.fD_type = fData::get_FDATA_TYPE( options.at(1) );
        // Obtain the f data format.
        tarFData.fD_format = fData::get_FDATA_FORMAT( options.at(2) );
        // Obtain the input impedance of the measurement.
        if( options.at(3) == "R" ){
            tarFData.systImp = std::stod( options.at(4) );
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
        tarFData.Xr_vec.reInit( tarFData.IOcnt[0], tarFData.IOcnt[1], res_blk_size );
        tarFData.Xi_vec.reInit( tarFData.IOcnt[0], tarFData.IOcnt[1], res_blk_size );

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
    
            
            for( int i = 0; i < tarFData.IOcnt[0]; i++ ){
                for( int j = 0; j < tarFData.IOcnt[1]; j++ ){
                    // Read the next data mag.
                    iss >> word;    tmp_val = std::stod( word );
                    tarFData.Xr_vec.set( i, j, line_idx, tmp_val );
                    // Read the next data phase.
                    iss >> word;    tmp_val = std::stod( word );
                    tarFData.Xi_vec.set( i, j, line_idx, tmp_val );
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
    
        return;
    
}