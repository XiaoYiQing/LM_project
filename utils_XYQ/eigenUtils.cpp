#include "eigenUtils.h"




// ====================================================================== >>>>>
//      Write to File
// ====================================================================== >>>>>

void utils::VectorXd_to_file( const string& fileDir, const string& fileStem, 
    const Eigen::VectorXd& tarVec, int options ){

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
    unsigned int data_cnt = tarVec.size();

    // Initialize temporary complex variable.
    double tmp_cplx = 0;

    // Set the precision of the number being printed.
    unsigned int precision = 10;
    
// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Write Info Lines
// ---------------------------------------------------------------------- >>>>>

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Write Data Lines
// ---------------------------------------------------------------------- >>>>>

    for( unsigned int z = 0; z < data_cnt; z++ ){

        if( tarVec(z) >= 0 ){
            file << "+";
        }
        file << std::scientific << std::setprecision(precision) << tarVec(z);

        if( z < data_cnt - 1 ){
            file << "\n";
        }

    }

// ---------------------------------------------------------------------- <<<<<

    // Close the file once all is done.
    file.close();

}



void utils::Vd_to_file( const string& fileDir, const string& fileStem, 
        vector<double> tarVec, int options )
{

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
    unsigned int data_cnt = tarVec.size();

    // Initialize temporary complex variable.
    double tmp_cplx = 0;

    // Set the precision of the number being printed.
    unsigned int precision = 10;
    
// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Write Data Lines
// ---------------------------------------------------------------------- >>>>>

    for( unsigned int z = 0; z < data_cnt; z++ ){

        if( tarVec.at(z) >= 0 ){
            file << "+";
        }
        file << std::scientific << std::setprecision(precision) << tarVec.at(z);

        if( z < data_cnt - 1 ){
            file << "\n";
        }

    }

// ---------------------------------------------------------------------- <<<<<

    // Close the file once all is done.
    file.close();

}




vector<double> utils::file_to_Vd( const string& fullFileName ){

// ---------------------------------------------------------------------- >>>>>
//      Full Name Parsing
// ---------------------------------------------------------------------- >>>>>

    std::filesystem::path fullFilePath( fullFileName );
    std::string fileDir = fullFilePath.parent_path().string();
    std::string fileStem = fullFilePath.stem().string();
    std::string fileExt = fullFilePath.extension().string();
    
    if( fileExt != ".txt" && fileExt != ".csv" ){
        throw std::invalid_argument( "Target data file must have extension '.txt' or '.csv'." );
    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Data Read Prep
// ---------------------------------------------------------------------- >>>>>

    // Open the input file stream.
    std::ifstream inputFile( fullFilePath );
    if( !inputFile ){
        throw std::invalid_argument( "Failed to open stream on target data file." );
    }else{
        cout << "File \"" + fileStem + fileExt + "\": stream opened successful." << endl;
    }

    // The variables holding the line and word currently read, respectively.
    string line, word;
    // Initialize index variable for the row index.
    unsigned int row_idx = 0;
    // Boolean flag for indicating whether the line stream has reached the data lines.
    bool end_reached = false;
    // Initialize temporary double variables for storing values.
    double tmp_val = 0;

    // Initialize memory allocation control variables.
    unsigned int res_blk_size = 200;
    unsigned int curr_vec_size = 0;
    unsigned int curr_alloc_size = res_blk_size;
    // Initialize vector.
    vector<double> vec = vector<double>();
    vec.reserve( res_blk_size );

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Data Read
// ---------------------------------------------------------------------- >>>>>

    // Parse through the file lines.
    while( getline( inputFile, line ) ){

        // Update the current vector size.
        curr_vec_size++;

        // Set the stream for the current line.
        istringstream iss(line);
        // Read the first word.
        if( iss >> word ){

            // Translate the word into a double value freq.
            try{
                tmp_val = std::stod( word );
            }catch( ... ){
                inputFile.close();
                throw;
            }
            // Push the current value at the back of the vector.
            vec.push_back( tmp_val );

        }else{
            throw std::runtime_error( "Real vector data file must not have any empty line." );
        }

        // Second word read on the line, which is expected to fail.
        if( iss >> word ){
            throw std::runtime_error( "Real vector data file must not have more than 1 word per line." );
        }

        // Allocate more memory if need be.
        if( curr_vec_size >= curr_alloc_size ){

            vec.reserve( curr_vec_size + res_blk_size );
            curr_alloc_size += res_blk_size;

        }

    }

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Function End Cleanup
// ---------------------------------------------------------------------- >>>>>

    // End file stream
    inputFile.close();
    // Downsize the vector.
    vec.shrink_to_fit();

    // End of function messages.
    cout << "File \"" + fileStem + fileExt + "\": successfully read." << endl;
    cout << "Vector size: " << curr_vec_size << endl;
    // cout << endl;

// ---------------------------------------------------------------------- <<<<<

    return vec;

}


void utils::MatrixXd_to_file( const string& fileDir, const string& fileStem, 
    const Eigen::MatrixXd& tarMat, int options )
{

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
    unsigned int row_cnt = tarMat.rows();
    unsigned int col_cnt = tarMat.cols();

    // Initialize temporary complex variable.
    double tmp_cplx = 0;

    // Set the precision of the number being printed.
    unsigned int precision = 10;
    
// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Write Data Lines
// ---------------------------------------------------------------------- >>>>>

    for( unsigned int i = 0; i < row_cnt; i++ ){
        for( unsigned int j = 0; j < col_cnt; j++ ){

            if( tarMat( i, j ) >= 0 ){
                file << "+";
            }
            file << std::scientific << std::setprecision(precision) << tarMat( i, j );

            if( j < col_cnt - 1 ){
                file << ", ";
            }

        }

        if( i < row_cnt - 1 ){
            file << "\n";
        }

    }

// ---------------------------------------------------------------------- <<<<<

    // Close the file once all is done.
    file.close();

}

// ====================================================================== <<<<<
