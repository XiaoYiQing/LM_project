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

    /*
    Write a VectorXd variable's content into a text file.
    The double values are simply placed one after each other on consecutive rows.
    Each row has only one value printed.
    */
    void VectorXd_to_file( const string& fileDir, const string& fileStem, 
        const Eigen::VectorXd& tarVec, int options );
    
    /*
    Write a vector<double> variable's content into a text file.
    The double values are simply placed one after each other on consecutive rows.
    Each row has only one value printed.
    */ 
    void Vd_to_file( const string& fileDir, const string& fileStem, 
        vector<double> tarVec, int options );

    /*
    Read an eigen vector saved in a .txt or .csv file and translate it into a 
    eigen vector variable.
    */
    vector<double> file_to_Vd( const string& fullFileName );

    // TODO: real mat to file.
    void MatrixXd_to_file( const string& fileDir, const string& fileStem, 
        const Eigen::MatrixXd& tarVec, int options );

    // TODO: file to real mat.
    Eigen::MatrixXd file_to_MatrixXd( const string& fullFileName );

// ====================================================================== <<<<<

}




#endif  // EIGENUTILS