#include "eigenUtils.h"





Eigen::MatrixXd utils::gen_rand_MatrixXd( unsigned int row_cnt, unsigned int col_cnt ){

    // Initialize random matrix.
    Eigen::MatrixXd randMat( row_cnt, col_cnt );

    for( unsigned int i = 0; i < row_cnt; i++ ){
        // Generate the current row.
        vector<double> curr_row = *utils::rDoubleGen( -1, 1, col_cnt );
        // Insert the current row.
        for( unsigned int j = 0; j < col_cnt; j++ ){
            randMat(i,j) = curr_row[j];
        }
    }

    return randMat;

}

// ====================================================================== >>>>>
//      Write to File
// ====================================================================== >>>>>

void utils::VectorXd_to_file( const string& fileDir, const string& fileStem, 
    const Eigen::VectorXd& tarVec, unsigned int decimCnt ){

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
    unsigned int precision = decimCnt;
    
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
        vector<double> tarVec, unsigned int decimCnt )
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
    unsigned int precision = decimCnt;
    
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
    const Eigen::MatrixXd& tarMat, unsigned int decimCnt )
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
    unsigned int precision = decimCnt;
    
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
                file << " ";
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




Eigen::MatrixXd utils::file_to_MatrixXd( const string& fullFileName ){

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
    unsigned int col_idx = 0;
    // Boolean flag for indicating whether the line stream has reached the data lines.
    bool end_reached = false;
    // Initialize temporary double variables for storing values.
    double tmp_val = 0;

    // Initialize memory allocation control variables.
    unsigned int res_blk_size = 200;
    unsigned int curr_vec_size = 0;
    unsigned int curr_alloc_size = res_blk_size;

    

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Number of Lines and Columns Count
// ---------------------------------------------------------------------- >>>>>

    // Initialize row and col counts.
    unsigned int row_cnt = 0;
    unsigned int col_cnt = 0;

    // Obtain the first line.
    getline( inputFile, line ); row_cnt++;

    // Define regex with space or comma as delimiter
    std::regex delimiter("[ ,]+");
    // Split the first line with space and/or commas.
    std::sregex_token_iterator iter(line.begin(), line.end(), delimiter, -1);
    std::sregex_token_iterator end;

    vector<string> lineSplit;
    while (iter != end) {
        lineSplit.push_back(*iter++);
    }
    col_cnt = lineSplit.size();

    // // Set the stream for the current line.
    // istringstream iss(line);
    // while( iss >> word ){
    //     col_cnt++;
    // }

    // First pass: count lines
    while ( getline( inputFile, line ) ) {
        row_cnt++;
    }

    // Move back to beginning
    inputFile.clear();                  // clear EOF flag
    inputFile.seekg(0, std::ios::beg);  // go back to start

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Data Read
// ---------------------------------------------------------------------- >>>>>

    // Initialize the data matrix.
    Eigen::MatrixXd datMat( row_cnt, col_cnt );

    // Parse through the file lines.
    for( unsigned int i = 0; i < row_cnt; i++ ){

        // Obtain the current line.
        getline( inputFile, line );

        // Split the first line with space and/or commas.
        std::sregex_token_iterator iter(line.begin(), line.end(), delimiter, -1);
        std::sregex_token_iterator end;

        unsigned int j = 0;
        while (iter != end) {

            if( j >= col_cnt ){
                throw runtime_error( "file_to_MatrixXd: inconsistent number of columns across the rows in the data file." );
            }
            // Try parsing the next word into double and store into data matrix.
            try{
                datMat(i,j) = std::stod(*iter++);
            }catch(...){
                inputFile.close();
                throw;
            }

            j++;

        }
        

        // // Set the stream for the current line.
        // istringstream iss(line);

        // for( unsigned int j = 0; j < col_cnt; j++ ){

        //     try{
                
        //         // Read the next word.
        //         iss >> word;
        //         // Translate the word into a double value freq.
        //         readMat(i,j) = std::stod( word );

        //     }catch( ... ){
                
        //         inputFile.close();
        //         throw;

        //     }

        // }
        
    }

// ---------------------------------------------------------------------- <<<<<

    return datMat;

}



// ====================================================================== <<<<<
