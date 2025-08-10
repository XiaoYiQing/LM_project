#include "eigenUtils.h"




// ====================================================================== >>>>>
//      Write to File
// ====================================================================== >>>>>

void utils::vec_to_file( const string& fileDir, const string& fileStem, 
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

// ====================================================================== <<<<<
