#ifndef EIGENUTILS_H
#define EIGENUTILS_H

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <string>




using namespace std;



namespace utils{

// ====================================================================== >>>>>
//      Write to File
// ====================================================================== >>>>>

    void vec_to_file( const string& fileDir, const string& fileStem, 
        const Eigen::VectorXd& tarVec, int options );

    /*
    Read an eigen vector saved in a .txt or .csv file and translate it into a 
    eigen vector variable.
    */
    vector<double> file_to_vec( const string& fullFileName );

    // TODO: real mat to file.

    // TODO: file to real mat.

// ====================================================================== <<<<<

}




#endif  // EIGENUTILS