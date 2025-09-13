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
        const Eigen::VectorXd& tarVec, unsigned int decimCnt );
    
    /*
    Write a vector<double> variable's content into a text file.
    The double values are simply placed one after each other on consecutive rows.
    Each row has only one value printed.
    */ 
    void Vd_to_file( const string& fileDir, const string& fileStem, 
        vector<double> tarVec, unsigned int decimCnt );

    /*
    Read an vector<double> saved in a .txt or .csv file and translate it into a 
    eigen vector variable.
    NOTE: expect simplest data format in the file, which is simply entry of each
    row placed on their own line consecutively.
    */
    vector<double> file_to_Vd( const string& fullFileName );

    /*
    Write a Eigen::MatrixXd variable's content into a text file.
    Entries of the matrix are written into the file preserving a semblant matrix 
    format. Columns are separated by single spaces and rows are separated by lines.
    */
    void MatrixXd_to_file( const string& fileDir, const string& fileStem, 
        const Eigen::MatrixXd& tarVec, unsigned int decimCnt );

    /*
    Read an Eigen::MatrixXd saved in a .txt or .csv file and translate it into a 
    eigen vector variable.
    */
    Eigen::MatrixXd file_to_MatrixXd( const string& fullFileName );

// ====================================================================== <<<<<

}




#endif  // EIGENUTILS