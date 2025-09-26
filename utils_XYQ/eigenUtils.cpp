#include "eigenUtils.h"






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


void utils::Vd_to_file_sci( const string& fileDir, const string& fileStem, 
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



void utils::MatrixXcd_to_file( const string& fileDir, const string& fileStem, 
    const Eigen::MatrixXcd& tarMat, unsigned int decimCnt )
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

            if( tarMat( i, j ).real() >= 0 ){
                file << "+";
            }
            file << std::scientific << std::setprecision(precision) << tarMat( i, j ).real();
            
            file << " ";

            if( tarMat( i, j ).imag() >= 0 ){
                file << "+";
            }
            file << std::scientific << std::setprecision(precision) << tarMat( i, j ).imag();

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


Eigen::MatrixXcd utils::file_to_MatrixXcd( const string& fullFileName ){

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

    // Initialize temporary double variables for storing values.
    double real_ij = 0;
    double imag_ij = 0;

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Number of Lines and Columns Count
// ---------------------------------------------------------------------- >>>>>

    // Initialize row and col counts for the data file.
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
    if( col_cnt % 2 != 0 ){
        throw runtime_error( "file_to_MatrixXcd: Number of data columns is not even." );
    }

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
    
    // Define number of data matrix rows and columns.
    unsigned int data_row_cnt = row_cnt;
    unsigned int data_col_cnt = col_cnt/2;

    // Initialize the data matrix.
    Eigen::MatrixXcd datMat( data_row_cnt, data_col_cnt );

    // Parse through the file lines.
    for( unsigned int i = 0; i < data_row_cnt; i++ ){

        // Obtain the current line.
        getline( inputFile, line );

        // Split the first line with space and/or commas.
        std::sregex_token_iterator iter(line.begin(), line.end(), delimiter, -1);
        std::sregex_token_iterator end;

        unsigned int j = 0;
        while (iter != end) {

            if( j >= data_col_cnt ){
                throw runtime_error( "file_to_MatrixXd: inconsistent number of columns across the rows in the data file." );
            }
            // Try parsing the next word into double and store into data matrix.
            try{
                real_ij = std::stod(*iter++);
                imag_ij = std::stod(*iter++);
                datMat(i,j) = complex<double>( real_ij, imag_ij );
            }catch(...){
                inputFile.close();
                throw;
            }

            j++;

        }
        
    }

// ---------------------------------------------------------------------- <<<<<

    return datMat;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Support Functions
// ====================================================================== >>>>>


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


Eigen::MatrixXd utils::gen_randn_MatrixXd( unsigned int row_cnt, unsigned int col_cnt ){

    std::random_device rd;  // Seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> dist(0, 1);

    Eigen::MatrixXd newMat( row_cnt, col_cnt );

    for ( unsigned int i = 0; i < row_cnt; ++i ){
        for ( unsigned int j = 0; j < col_cnt; ++j ){
            newMat(i, j) = dist(gen);
        }
    }

    return newMat;

}


Eigen::MatrixXd utils::gen_orth_basis( Eigen::MatrixXd tarMat ){

    // Generate QR-decomposition out of the target matrix.
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(tarMat);

    // The columns of Q constitute an orthonormal basis.
    Eigen::MatrixXd Q = qr.householderQ();

    return Q;

}


utils::SVD_econ::SVD_econ(){

    this->U = Eigen::MatrixXd();
    this->S = Eigen::MatrixXd();
    this->V = Eigen::MatrixXd();

}

utils::SVD_econ::SVD_econ( Eigen::MatrixXd& tarMat ){
    
    unsigned int k = min( tarMat.rows(), tarMat.cols() );

    Eigen::JacobiSVD<Eigen::MatrixXd> svd( tarMat, Eigen::ComputeThinU | Eigen::ComputeThinV );

    // Extract truncated components
    this->U = svd.matrixU().leftCols(k);
    this->S = svd.singularValues().head(k);
    this->V = svd.matrixV().leftCols(k);

}

utils::rSVD::rSVD(){
    Uk = Eigen::MatrixXd();
    Sk = Eigen::VectorXd();
    Vk = Eigen::MatrixXd();
}

utils::rSVD::rSVD( Eigen::MatrixXd& A, unsigned int k ){

    // Obtain the number of rows and columns of the target matrix.
    unsigned int A_rows = A.rows();
    unsigned int A_cols = A.cols();

    /*
    Select the size of the solution space. 
    Typically, it is chosen as twice as the required number k of singular values
    in order to make sure the k first singular values are properly approximated.
    Of course, this size cannot be more than the total available number of singular values.
    */
    unsigned int P = min( 2*k, A_cols );

    /*
    Generate a normally distributed random matrix,
    imprint it on the target matrix,
    generate orthonormal basis from this solution space.
    */
    Eigen::MatrixXd X = utils::gen_randn_MatrixXd( A_cols, P );
    Eigen::MatrixXd Y = A*X;
    Eigen::MatrixXd W1 = utils::gen_orth_basis( Y );
    
    Eigen::MatrixXd B = W1.adjoint() * A;

    // Generate the approximate SVD solution.
    utils::SVD_econ svd_econ_res ( B );
    Eigen::MatrixXd W2 = svd_econ_res.U;
    Eigen::MatrixXd U = W1 * W2;

    // Update k to a number within limit, if need be.
    unsigned int tmp = U.cols();
    k = min( k, tmp );

    unsigned int U_row_cnt = U.rows();
    unsigned int V_row_cnt = svd_econ_res.V.rows();

    Uk = U.block( 0, 0, U_row_cnt, k );
    Sk = svd_econ_res.S.segment( 0, k );
    Vk = svd_econ_res.V.block( 0, 0, V_row_cnt, k );

}

// ====================================================================== <<<<<
