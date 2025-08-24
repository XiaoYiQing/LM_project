#ifndef FDATA_H
#define FDATA_H

#include <cmath>
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <magic_enum.hpp>
#include <regex>
#include <sstream>
#include <string>
#include <vector>  

#include "Matrix3DXd.h"
#include "numUtils.h"

using namespace std;


class fData{    

public:

// ====================================================================== >>>>>
//      Data From File Utility
// ====================================================================== >>>>>

    /*
    Function to retrieve frequency data from files ending with the extension of format
    ".sXp" where X is any positive integer representing number of I/O ports.
    */
    static void read_sXp_file( fData& tarFData, const string& fullFileName );

    /*
    Read the specific 2-port S-parameter data file, which follows a different parsing
    rule than 3 or higher number of ports S data files.
    */
    static void read_s2p_file( fData& tarFData, const string& fullFileName );


    /*
    Function to retrieve freq. data from S-parameter file data generated from the 
    LTspice software.
    NOTE:
    - LTspice only provides two-port network S-parameters by default (That is S11, S12, 
        S21, S22), so the parser will only extract 2x2 S-parameter matrices.
    - This parser only reads data given in the (dB,degree (non-radian)) format.
    */
    static void read_LTspice_Sp_file( fData& tarFData, const string& fullFileName );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Class Enum "METRIC_PREFIX" Help Functions
// ====================================================================== >>>>>

    /*
    Enum representing the metric system prefix symbols.
    In the order defined:
        pico, nano, micro, micro (alt), milli, centi, deci, deca, hecto, kilo, mega, giga, tera
    */
    enum class METRIC_PREFIX{ NONE, p, n, μ, mu, m, c, d, da, h, k, M, G, T };

    // The number of enum entries in the enum "METRIC_PREFIX" (Uses magic enum).
    const static int METRIC_PREFIX_Count = (int) magic_enum::enum_count<METRIC_PREFIX>();

    // Obtain the string of the target enum case (Uses magic enum).
    static string get_METRIC_PREFIX_Str( METRIC_PREFIX tar_METRIC_PREFIX );
    // Obtain the enum matching the enum integer index.
    static METRIC_PREFIX get_METRIC_PREFIX_AtIdx( int idx );
    // Obtain the prefix using the equivalent string symbol.
    static METRIC_PREFIX get_METRIC_PREFIX( string strSymbol );
    // Obtain the numerical value of the metric prefix.
    static double get_METRIC_PREFIX_val( METRIC_PREFIX tar_METRIC_PREFIX );
    // Obtain the next higher prefix TODO
    static METRIC_PREFIX get_higher_prefix( METRIC_PREFIX tar_METRIC_PREFIX );

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

    /*
    Initialization with known frequency vector and frequency data.
    NOTE: scale, data-type, and data format are all set to NONE by default. You 
    have to change them by yourself afterward initialization.
    */
    fData( Eigen::VectorXd&, Matrix3DXd&, Matrix3DXd& );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Data Editing
// ====================================================================== >>>>>

    /*
    Re-initialize the fData by setting the f vector and the f data vector to
    only zeros values with the given number of ports and frequency points.
    This function does not affect format settings.
    */
    void reInit( unsigned int out_cnt, unsigned int in_cnt, unsigned int f_cnt );

    /*
    Set the data of the target objects to match those of the reference object.
    NOTE: Data copying does NOT perform settings copying such as frequency data type and format.
    */
    static void copy_data( fData& tarObj, const fData& refObj );

    /*
    Set the settings of the target objects to match those of the reference object.
    NOTE: Settings copying does NOT perform any kind of data format conversion.
    The target will adopt all settings of the reference directly without any check
    and modifition on the data, so the frequency data and the frequency array are 
    not touched.
    */
    static void copy_settings( fData& tarObj, const fData& refObj );

    /*
    Switch the format of the data from the current format to the specified new format.
    */
    void data_format_Switch( FDATA_FORMAT newFormat );
    /*
    Switch the data prefix, which applies the corresponding rescaling to the f vector.
    */
    void data_prefix_switch( METRIC_PREFIX newPref );

    /*
    Normalize the frequency array. NOTE: this function simply selects an appropriate
    metric prefix switch such TODO.
    */
    void f_normalize();

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Support Functions
// ====================================================================== >>>>>

    /*
    Return true if the zero frequency point is within the f vector.
    NOTE: will only check the first entry of the f vector, which is the only
    slot where the DC point is allowed to be.
    */
    bool hasDC() const;

    /*
    Create a reduced set of frequency data using a smaller linearly distributed
    index set.
    */
    shared_ptr<fData> red_partit_lin( unsigned int rSize );
    
    /*
    Create a reduced set of frequency data using a subset index array.
    */
    shared_ptr<fData> red_partit( vector< unsigned int > fr_idx_vec );

    /*
    Create two partitions from the original frequency data set.
    The partitions are decided by interleaving indexing.
    Partition 1 always has the first frequency entry.
    */
    vector< shared_ptr<fData> > gen_2_partit() const;

    /*
    Create two index arrays having mutually exclusive indices which
    can serve to create two partitions from the existing f data set.
    */
    vector< vector< unsigned int > > gen_2_partit_idx_arr() const;

    /*
    Create two partitions from the original frequency data set.
    Partition 1 is decided by the input index array.
    Partition 2 receives all leftover.
    TODO: verify if partition 1 index array have repeat.
    */
    vector< shared_ptr<fData> > gen_2_partit( const vector< unsigned int >& p1_idx );

    /*
    Generate a complex conjugate set of data.
    This function assumes the frequency data are ordered in the ascending order
    of frequency magnitudes (If DC point is present, it MUST be first entry).
    The generated freq. data set is a complement set, so it will not include
    the DC data point if it is present in the original.
    */
    fData gen_cplx_conj_set() const;

    /*
    Generate a copy of the present fData object but with complex conjugate frequency
    and data inserted in interleaving fashion.
    This means each frequency and the corresponding frequency data is immediately 
    followed by their complex conjugate in the arrays.
    */
    shared_ptr<fData> gen_cplx_conj_comb() const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>

    /*
    Set the number of inputs and outputs.
    NOTE: This function expunges existing frequency data.
    */
    void set_IO_cnt( unsigned int new_in_cnt, unsigned int new_out_cnt );

    /*
    Set the number of outputs.
    NOTE: This function expunges existing frequency data.
    */
    void set_out_cnt( unsigned int new_out_cnt );
    // Obtain the number of outputs (Number of rows in a data matrix).
    unsigned int get_out_cnt() const;
    /*
    Set the number of inputs.
    NOTE: This function expunges existing frequency data.
    */
    void set_in_cnt( unsigned int new_in_cnt );
    // Obtain the number of inputs (Number of columns in a data matrix).
    unsigned int get_in_cnt() const;
    // Obtain the number of frequency points.
    unsigned int get_f_cnt() const;
    // Obtain the string representation of the frequency scale or metric prefix.
    string get_f_scale_str() const;
    // Obtain the numerical value of the frequency scale or metric prefix.
    double get_f_scale_num() const;

    // Set the frequency value at the target frequency index.
    void set_fval_at( unsigned int f_idx, double f_val );
    // Obtain the frequency value at target frequency index.
    double get_fval_at( unsigned int f_idx ) const;

    /*
    Replace the block of continuous frequency within the frequency vector with the given block.
    - lead: the starting index of the block.
    - f_blk: the block of frequency to insert.
    */
    void set_fval_block( unsigned int lead, const Eigen::VectorXd& f_blk );

    // Obtain the complex frequency value at target frequency index.
    // TODO: Consider putting the 2*pi scaling in the future.
    std::complex<double> get_cplx_f_at( unsigned int f_idx ) const;

    // Set the real part data matrix at target frequency index.
    void set_reData_at_f( unsigned int f_idx, const Eigen::MatrixXd& new_rePart );
    // Obtain the real part data matrix at target frequency index.
    Eigen::MatrixXd get_reData_at_f( int f_idx ) const;

    // Set the real part data matrix at target frequency index.
    void set_imData_at_f( unsigned int f_idx, const Eigen::MatrixXd& new_imPart );
    // Obtain the imaginary part data matrix at target frequency index.
    Eigen::MatrixXd get_imData_at_f( unsigned int f_idx ) const;

    // Set the complex data at the target frequency index.
    void set_cplxData_at_f( unsigned int f_idx, Eigen::MatrixXcd& new_mat );
    // Obtain the complex data matrix at the target frequency index.
    Eigen::MatrixXcd get_cplxData_at_f( unsigned int f_idx ) const;

    /*
    Replace the block of continuous data within the frequency data array with the given block.
    - lead: the starting index of the block.
    - newBlk_re: insert block real part.
    - newBlk_im: insert block imag part.
    */
    void set_cplxData_block( unsigned int lead, const Matrix3DXd& newBlk_re, const Matrix3DXd& newBlk_im );

    // Obtain the f vector.
    Eigen::VectorXd getF_vec() const;
    // Obtain the frequency data part A (Real part or magnitude).
    Matrix3DXd getXr_vec() const;
    // Obtain the frequency data part B (Imaginary part or phase).
    Matrix3DXd getXi_vec() const;



// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Printing
// ====================================================================== >>>>>

    void print_to( const string& fileDir, const string& fileStem, int options );

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
    Eigen::VectorXd f_vec;

    // The vector of frequency data matrices real part.
    // vector< Eigen::MatrixXd > Xr_vec;
    Matrix3DXd Xr_vec;
    // The vector of frequency data matrices imaginary part.
    // vector< Eigen::MatrixXd > Xi_vec;
    Matrix3DXd Xi_vec;


private:


};






#endif  // FDATA_H




