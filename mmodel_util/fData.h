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
    In the order defined:
        pico, nano, micro, milli, centi, deci, deca, hecto, kilo, mega, giga, tera
    */
    enum class METRIC_PREFIX{ NONE, p, n, μ, m, c, d, da, h, k, M, G, T };

    // The number of enum entries in the enum "METRIC_PREFIX" (Uses magic enum).
    const static int METRIC_PREFIX_Count = (int) magic_enum::enum_count<METRIC_PREFIX>();

    // Obtain the string of the target enum case (Uses magic enum).
    static string get_METRIC_PREFIX_Str( METRIC_PREFIX tar_METRIC_PREFIX );
    // Obtain the enum matching the enum integer index.
    static METRIC_PREFIX get_METRIC_PREFIX_AtIdx( int idx );
    // Obtain the prefix using the equivalent string symbol.
    static METRIC_PREFIX get_METRIC_PREFIX( string strSymbol );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Class Enum "FDATA_TYPE" Help Functions
// ====================================================================== >>>>>

    /*
    Enum representing the frequency data types.
    S stands for scattering paramreter, Y for admittance parameter, etc.
    */
    enum class FDATA_TYPE{ NONE, S, Y, Z, H, G };

    // The number of enum entries in the enum "FDATA_TYPE" (Uses magic enum).
    const static int FDATA_TYPE_Count = (int) magic_enum::enum_count<FDATA_TYPE>();

    // Obtain the string of the target enum case (Uses magic enum).
    static string get_FDATA_TYPE_Str( FDATA_TYPE tar_FDATA_TYPE );
    // Obtain the enum matching the enum integer index.
    static FDATA_TYPE get_FDATA_TYPE_AtIdx( int idx );
    // Obtain the prefix using the equivalent string symbol.
    static FDATA_TYPE get_FDATA_TYPE( string strSymbol );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Class Enum "FDATA_FORMAT" Help Functions
// ====================================================================== >>>>>
    
    /*
    Enum representing the frequency data types.
    In the order defined:
        DB: Decibel Magnitude/Degree Angle
        MA: Linear Magnitude/Degree Angle (polar form)
        RI: Real Part/Imaginary Part (rectangular form)
    */
    enum class FDATA_FORMAT{ NONE, DB, MA, RI };

    // The number of enum entries in the enum "FDATA_FORMAT" (Uses magic enum).
    const static int FDATA_FORMAT_Count = (int) magic_enum::enum_count<FDATA_FORMAT>();

    // Obtain the string of the target enum case (Uses magic enum).
    static string get_FDATA_FORMAT_Str( FDATA_FORMAT tar_FDATA_FORMAT );
    // Obtain the enum matching the enum integer index.
    static FDATA_FORMAT get_FDATA_FORMAT_AtIdx( int idx );
    // Obtain the prefix using the equivalent string symbol.
    static FDATA_FORMAT get_FDATA_FORMAT( string strSymbol );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

    fData();

// ====================================================================== <<<<<

protected:


    // Number of inputs and outputs, respectively.
    int IOcnt[2];

    // The metric prefix for the frequency vector.
    fData::METRIC_PREFIX f_pref;
    // The frequency data type, such as S-parameter, Y-parameter, etc.
    fData::FDATA_TYPE fD_type;
    // The frequency data format (linear mag/phase, log mag/phase, real/img).
    fData::FDATA_FORMAT fD_format;
    // The system impedence at which the measurements were taken.
    double systImp;

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




