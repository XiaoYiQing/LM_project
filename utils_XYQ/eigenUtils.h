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

    // TODO: file to vector.

    // TODO: real mat to file.

    // TODO: file to real mat.

// ====================================================================== <<<<<

}




#endif  // EIGENUTILS