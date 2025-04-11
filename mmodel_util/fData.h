#ifndef FDATA_H
#define FDATA_H


#include <Eigen/Dense>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <magic_enum.hpp>
#include <regex>
#include <sstream>
#include <string>   
#include <vector>  

using namespace std;

class fData{    

public:

    /*
    Function to retrieve frequency data from files ending with the extension of format
    ".sXp" where X is any positive integer representing number of I/O ports.
    */
    static void read_sXp_file( fData& tarFData, const string& fullFileName );


// ====================================================================== >>>>>
//      Class Enum "METRIC_PREFIX" Help Functions
// ====================================================================== >>>>>

    /*
    Enum representing the metric system prefix symbols.
    */
    enum class METRIC_PREFIX{ p, n, mu, m, c, d, da, h, k, M, G, T };

    // The number of enum entries in the enum "CHK_PIECE" (Uses magic enum).
    const static int METRIC_PREFIX_Count = (int) magic_enum::enum_count<METRIC_PREFIX>();

    // Obtain the string of the target enum case (Uses magic enum).
    static string get_METRIC_PREFIX_Str( METRIC_PREFIX tar_METRIC_PREFIX );
    // Obtain the enum matching the enum integer index.
    static METRIC_PREFIX get_METRIC_PREFIX_AtIdx( int idx );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

    fData();

// ====================================================================== <<<<<

protected:


    // Number of inputs and outputs, respectively.
    int IOcnt[2];

    /*
    The vector of frequencies in hertz (unit saved separately).
    */
    Eigen::VectorXd fr_vec;

    // The vector of frequency data matrices real part.
    vector< Eigen::MatrixXd > Xr_vec;
    // The vector of frequency data imaginary real part.
    vector< Eigen::MatrixXd > Xi_vec;

private:


};






#endif  // FDATA_H




