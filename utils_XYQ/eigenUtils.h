#ifndef EIGENUTILS_H
#define EIGENUTILS_H

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <string>

#include "numUtils.h"


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
        const Eigen::MatrixXd& tarMat, unsigned int decimCnt );

    /*
    Read an Eigen::MatrixXd saved in a .txt or .csv file and translate it into a 
    eigen vector variable.
    */
    Eigen::MatrixXd file_to_MatrixXd( const string& fullFileName );
    
    /*
    Write a Eigen::MatrixXd variable's content into a text file.
    Entries of the matrix are written into the file preserving a semblant matrix 
    format. Columns are separated by single spaces and rows are separated by lines.
    Each data has its real part placed on a former column and its imaginary part placed
    on the immediate following column, which means a data matrix of X columns generates
    a text data file of 2*X columns.
    */
    void MatrixXcd_to_file( const string& fileDir, const string& fileStem, 
        const Eigen::MatrixXcd& tarMat, unsigned int decimCnt );

    /*
    Read an Eigen::MatrixXcd saved in a .txt or .csv file and translate it into a 
    eigen vector variable.
    */
    Eigen::MatrixXcd file_to_MatrixXcd( const string& fullFileName );


// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Support Functions
// ====================================================================== >>>>>

    /*
    Generate random "row_cnt" by "col_cnt" MatrixXd with entries within the range of [-1,+1].
    */
    Eigen::MatrixXd gen_rand_MatrixXd( unsigned int row_cnt, unsigned int col_cnt );

    /*
    Function generates an orthonormal basis for the target matrix.
    */
    Eigen::MatrixXd gen_orth_basis( Eigen::MatrixXd tarMat );

// ====================================================================== <<<<<

}




#endif  // EIGENUTILS