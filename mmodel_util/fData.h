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

/**
 * Template object for storing frequency data.
 * 
 * @note 
 * 1- DC point is ALWAYS at index 0.
 * 2- No repeated frequency entries.
 */
class fData{    

public:

// ====================================================================== >>>>>
//      Data From File Utility
// ====================================================================== >>>>>

    /**
     * Function to retrieve frequency data from files ending with the extension of format
     * ".sXp" where X is any positive integer representing number of I/O ports.
     * 
     * @param tarFData The fData object into which the read data is inserted.
     * @param fullFileName The complete file name (directory, file stem, file extension) of 
     * the .sXp file to be read.
     */
    static void read_sXp_file( fData& tarFData, const string& fullFileName );

    /**
     * Read the specific 2-port S-parameter data file, which follows a different parsing
     * rule than 3 or higher number of ports S data files.
     * 
     * @param tarFData The fData object into which the read data is inserted.
     * @param fullFileName The complete file name (directory, file stem, file extension) of 
     * the .s2p file to be read.
     */
    static void read_s2p_file( fData& tarFData, const string& fullFileName );


    /**
     * Function to retrieve freq. data from S-parameter file data generated from the 
     * LTspice software.
     * 
     * @param tarFData The fData object into which the read data is inserted.
     * @param fullFileName The complete file name (directory, file stem, file extension) of 
     * the .s2p file to be read.
     * 
     * @NOTE:
     * 
     * - LTspice only provides two-port network S-parameters by default (That is S11, S12, 
     *   S21, S22), so the parser will only extract 2x2 S-parameter matrices.
     * 
     * - This parser only reads data given in the (dB,degree (non-radian)) format.
     * 
     * - LTspice does not print additional information such as input impedance. As such,
     * such values are set to default (Rin = 50, for example).
    */
    static void read_LTspice_Sp_file( fData& tarFData, const string& fullFileName );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Class Enum "METRIC_PREFIX" Help Functions
// ====================================================================== >>>>>

    /**
     * Enum representing the metric system prefix symbols.
     * In the order defined:
     * 
     * - pico, nano, micro, micro (alt), milli, centi, deci, deca, hecto, kilo, mega, giga, tera
     */
    enum class METRIC_PREFIX{ p, n, μ, mu, m, c, d, NONE, da, h, k, M, G, T };

    // The number of enum entries in the enum "METRIC_PREFIX" (Uses magic enum).
    const static int METRIC_PREFIX_Count = (int) magic_enum::enum_count<METRIC_PREFIX>();

    /**
     * Obtain the string of the target enum case (Uses magic enum).
     * 
     * @param tar_METRIC_PREFIX The target metric prefix.
     * @return The string representation of the target metric prefix.
     */
    static string get_METRIC_PREFIX_Str( METRIC_PREFIX tar_METRIC_PREFIX );
    /**
     * Obtain the enum matching the enum integer index (return -1 if failed).
     * 
     * @param idx Index of the target metric prefix.
     * @return The metric prefix associated to the target index.
     */
    static METRIC_PREFIX get_METRIC_PREFIX_AtIdx( int idx );
    /**
     * Obtain the prefix using the equivalent string symbol.
     * 
     * @param strSymbol The target string representation of a metric prefix.
     * @return The actual metric prefix represented by the target string.
     */
    static METRIC_PREFIX get_METRIC_PREFIX( string strSymbol );
    /**
     * Obtain the numerical value of the metric prefix.
     * 
     * @param tar_METRIC_PREFIX The target metric prefix.
     * @return The numerical value equivalent of one unit of the target metric prefix.
     */
    static double get_METRIC_PREFIX_val( METRIC_PREFIX tar_METRIC_PREFIX );
    
    /**
     * Obtain the next prefix:
     * - If "higher": next higher prefix.
     * - If not "higher": next lower prefix.
     * 
     * @param tar_METRIC_PREFIX The target metrix prefix.
     * @param higher The boolean indicating whether to go higher than the target prefix.
     * @return The next metrix prefix following the target one.
     */
    static METRIC_PREFIX get_METRIC_PREFIX_next( METRIC_PREFIX tar_METRIC_PREFIX, bool higher );
    /**
     * Obtain the appropriate prefix for the given value.
     * For example, inputting 123456 would return "k" for kilo whereas 
     * inputting 1234567 would return "M" for mega.
     * 
     * @param tarVal The numerical value for which a metric prefix is sought for.
     * @return The metric prefix most suited for the given numerical value.
     */
    static METRIC_PREFIX get_METRIC_PREFIX_for_val( double tarVal );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Class Enum "FDATA_TYPE" Help Functions
// ====================================================================== >>>>>

    /**
     * Enum representing the frequency data types.
     * S stands for scattering paramreter, Y for admittance parameter, etc.
     */
    enum class FDATA_TYPE{ NONE, S, Y, Z, H, G };

    /**
     * The number of enum entries in the enum "FDATA_TYPE" (Uses magic enum).
     */
    const static int FDATA_TYPE_Count = (int) magic_enum::enum_count<FDATA_TYPE>();

    /**
     * Obtain the string of the target enum case (Uses magic enum).
     * 
     * @param tar_FDATA_TYPE The target data type.
     * @return The string representation of the target data type.
     */
    static string get_FDATA_TYPE_Str( FDATA_TYPE tar_FDATA_TYPE );
    /**
     * Obtain the enum matching the enum integer index.
     * 
     * @param idx The target data type index.
     * @return The data type associated to the target index.
     */
    static FDATA_TYPE get_FDATA_TYPE_AtIdx( int idx );
    /**
     * Obtain the data type equivalent to the target string representation.
     * 
     * @param strSymbol The target data type string representation.
     * @return The data type associated to the string representation.
     */
    static FDATA_TYPE get_FDATA_TYPE( string strSymbol );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Class Enum "FDATA_FORMAT" Help Functions
// ====================================================================== >>>>>
    

    /**
     * Enum representing the frequency data types.
     * In the order defined:
     *  - DB: Decibel Magnitude/Degree Angle
     *  - MA: Linear Magnitude/Degree Angle (polar form)
     *  - RI: Real Part/Imaginary Part (rectangular form)
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
//      General Static Support Functions
// ====================================================================== >>>>>

    /**
     * Create two index arrays having mutually exclusive indices which 
     * can serve to create two partitions from the existing original set.
     */
    // static vector< vector< unsigned int > > gen_2_partit_idx_arr( unsigned int origSize );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

    fData();

    /**
     * Initialization with known frequency vector and frequency data.
     * 
     * @param f_vec Vector of frequency magnitudes (Hz).
     * @param Xr_vec Frequency data real part.
     * @param Xi_vec Frequency data imaginary part.
     * 
     * @NOTE: scale, data-type, and data format are all set to NONE by default. You 
     * have to change them by yourself afterward initialization.
     */
    fData( Eigen::VectorXd& f_vec, Matrix3DXd& Xr_vec, Matrix3DXd& Xi_vec );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Data Editing
// ====================================================================== >>>>>

    /**
     * Re-initialize the fData by setting the f vector and the f data vector 
     * to only zeros values with the given number of ports and frequency points. 
     * This function does not affect format settings.
     * 
     * @param out_cnt The number of outputs (Number of rows of data matrices).
     * @param in_cnt The number of inputs (Number of columns of data matrices).
     * @param f_cnt Number of frequency points (Number of data matrices).
     */
    void reInit( unsigned int out_cnt, unsigned int in_cnt, unsigned int f_cnt );

    /**
     * Set the data of the target objects to match those of the reference object.
     * This includes frequency array copying.
     * 
     * @param tarObj The fData instance where the data is to be copied into.
     * @param refObj The fData instance whose data is the source of the copying.
     * 
     * @note Data copying does NOT perform settings copying such as frequency data type and format.
     */
    static void copy_data( fData& tarObj, const fData& refObj );

    /**
     * Set the settings of the target objects to match those of the reference object.
     * 
     * @param tarObj The fData instance where the settings is to be copied into.
     * @param refObj The fData instance whose settings are the source of the copying.
     * 
     * @note Settings copying does NOT perform any kind of data format conversion. 
     * The target will adopt all settings of the reference directly without any check 
     * and modifition on the data, so the frequency data and the frequency array are 
     * not touched.
     */
    static void copy_settings( fData& tarObj, const fData& refObj );

    /**
     * Switch the format of the data from the current format to the specified new format.
     * 
     * @param newFormat The new data format to which the current data is to switch to.
     */
    void data_format_Switch( FDATA_FORMAT newFormat );

    /**
     * Switch the data prefix, which applies the corresponding rescaling to the f vector.
     * 
     * @param newPref The new metric prefix to which the current data is to rescale to.
     */
    void data_prefix_switch( METRIC_PREFIX newPref );

    /**
     * Normalize the frequency array with respect to known metric prefixes.
     * 
     * @note This function simply selects an appropriate metric prefix to switch to.
     */
    void f_rescale();

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Specialized Support Functions
// ====================================================================== >>>>>

    /**
     * Determine if the data has any 0 frequency (DC point) data point.
     * 
     * @return Boolean indicating whether the DC point is present in this instance of fData.
     * 
     * @note Will only check the first entry of the f vector, which is the only 
     * slot where the DC point is allowed to be.
     */
    bool hasDC() const;

    /**
     * Create a reduced set of frequency data using a smaller linearly distributed
     * index set.
     * 
     * @param rSize The size of the reduced linearly distributed frequency set.
     * @return The fData formed from the reduced set of frequency points.
     */
    fData red_partit_lin( unsigned int rSize );
    
    /**
     * Create a reduced set of frequency data using a subset index array.
     * 
     * @param fr_idx_vec The index array for containing the indices of frequency 
     *  data to extract to create the reduced fData.
     * @return The reduced fData containing the frequency data from the original 
     *  fData at the specified frequency points.
     */
    fData red_partit( const vector< unsigned int >& fr_idx_vec ) const;

    /**
     * Create two partitions from the original frequency data set. The 
     * partitions are decided by interleaving indexing. 
     * 
     * @return Vector of 2 shared_ptr leading to the 2 fData partitions.
     * 
     * @note Partition 1 always has the first (index  = 0) frequency entry.
     */
    vector< shared_ptr<fData> > gen_2_partit() const;

    /**
     * Create two index arrays following the interleaving indexing pattern and 
     * having mutually exclusive indices which can serve to create two partitions 
     * from the existing f data set.
     * 
     * @return Vector of 2 index vector arrays.
     */
    vector< vector< unsigned int > > gen_2_partit_idx_arr() const;

    /*
    Create two partitions from the original frequency data set.
    Partition 1 is decided by the input index array.
    Partition 2 receives all leftover.
    */
    // vector< shared_ptr<fData> > gen_2_partit( const vector< unsigned int >& p1_idx );

    /**
     * Generate a complex conjugate set of data.
     * This function assumes the frequency data are ordered in the ascending order
     * of frequency magnitudes (If DC point is present, it MUST be first entry).
     * The generated freq. data set is a complement set, so it will not include
     * the DC data point if it is present in the original.
     * 
     * @return Create a complete complement set of complex conjugate frequency data.
     */
    fData gen_cplx_conj_set() const;

    /**
     * Generate a copy of the present fData object but with complex conjugate frequency 
     * and data inserted in interleaving fashion.
     * This means each frequency and the corresponding frequency data is immediately 
     * followed by their complex conjugate in the arrays.
     * 
     * @return Create a fData which is the original fData instance injected with its complete
     *  complex conjugate set of data.
     */
    fData gen_cplx_conj_comb() const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Functions
// ====================================================================== >>>>>

    /**
     * Set the number of inputs and outputs.
     * 
     * @param new_in_cnt The new number of inputs.
     * @param new_out_cnt The new number of outputs.
     * 
     * @warning This function expunges existing frequency data.
     */
    void set_IO_cnt( unsigned int new_in_cnt, unsigned int new_out_cnt );

    /**
     * Set the number of outputs.
     * 
     * @param new_out_cnt The new number of outputs.
     * 
     * @warning This function expunges existing frequency data.
     */
    void set_out_cnt( unsigned int new_out_cnt );
    /**
     * Obtain the number of outputs (Number of rows in a data matrix).
     * 
     * @return The number of outputs.
     */
    unsigned int get_out_cnt() const;
    /**
     * Set the number of inputs.
     * 
     * @param new_in_cnt The new number of inputs.
     * 
     * @warning This function expunges existing frequency data.
     */
    void set_in_cnt( unsigned int new_in_cnt );
    /**
     * Obtain the number of inputs (Number of columns in a data matrix).
     * 
     * @return The number of inputs.
     */
    unsigned int get_in_cnt() const;
    /**
     * Obtain the number of frequency points.
     * 
     * @return The number of frequency points (number of data matrices).
     */
    unsigned int get_f_cnt() const;
    /**
     * Obtain the string representation of the frequency scale or metric prefix.
     * 
     * @return The string representation of the current adopted frequency metric prefix.
     */
    string get_f_scale_str() const;
    /**
     * Obtain the numerical value of the frequency scale or metric prefix.
     * 
     * @return The numerical value of the current adopted frequency metric prefix.
     */
    double get_f_scale_num() const;
    /**
     * Obtain the string representation of the frequency data type.
     * 
     * @return The string representation of the current adopted frequency data type.
     */
    string get_f_data_type_str() const;
    /**
     * Obtain the string representation of the frequency data format.
     * 
     * @return The string representation of the current adopted frequency data format.
     */
    string get_f_data_format_str() const;

    /**
     * Obtain the system impedance.
     * 
     * @return The system impedance associated to this frequency data set.
     */
    double get_systImp() const;
    /**
     * Set the system impedance.
     * 
     * @param The new system impedance to be associated to this frequency data set.
     */
    void set_systImp( const double in_imp );

    /**
     * Set the frequency value at the target frequency index.
     * 
     * @param f_idx The index at which the frequency value is to be changed.
     * @param f_val The new value to be assigned at the target frequency index.
     */
    void set_fval_at( unsigned int f_idx, double f_val );
    /**
     * Get the frequency value at the target frequency index.
     * 
     * @param f_idx The index at which the frequency value is to be fetched.
     * 
     * @return The frequency value at the target index.
     */
    double get_fval_at( unsigned int f_idx ) const;

    /**
     * Replace the block of continuous frequency within the frequency vector with the given block.
     * 
     * @param lead The starting index of the block.
     * @param f_blk the block of frequency to insert.
     */
    void set_fval_block( unsigned int lead, const Eigen::VectorXd& f_blk );

    // TODO: Consider putting the 2*pi scaling in the future.
    /**
     * Obtain the complex frequency value at target frequency index.
     * 
     * @param f_idx The index at which the frequency value is to be fetched.
     * 
     * @return The frequency value at the target index.
     */
    std::complex<double> get_cplx_f_at( unsigned int f_idx ) const;

    /**
     * Set the data matrix real part at the target frequency index.
     * 
     * @param f_idx The index at which the frequency data real part is to be changed.
     * @param new_rePart The new real part to be applied at the target index.
     */
    void set_reData_at_f( unsigned int f_idx, const Eigen::MatrixXd& new_rePart );
    /**
     * Obtain the real part data matrix at the target frequency index.
     * 
     * @param f_idx The index at which the frequency data real part is to be fetched.
     * @return The real part of the data matrix at the target index.
     */
    Eigen::MatrixXd get_reData_at_f( int f_idx ) const;

    /**
     * Set the data matrix imaginary part at the target frequency index.
     * 
     * @param f_idx The index at which the frequency data imaginary part is to be changed.
     * @param new_imPart The new imaginary part to be applied at the target index.
     */
    void set_imData_at_f( unsigned int f_idx, const Eigen::MatrixXd& new_imPart );
    /**
     * Obtain the imaginary part data matrix at the target frequency index.
     * 
     * @param f_idx The index at which the frequency data imaginary part is to be fetched.
     * @return The imaginary part of the data matrix at the target index.
     */
    Eigen::MatrixXd get_imData_at_f( unsigned int f_idx ) const;

    /**
     * Set the data matrix at the target frequency index.
     * 
     * @param f_idx The index at which the frequency data is to be changed.
     * @param new_mat The new data matrix to be applied at the target index.
     */
    void set_cplxData_at_f( unsigned int f_idx, Eigen::MatrixXcd& new_mat );
    /**
     * Obtain the data matrix at the target frequency index.
     * 
     * @param f_idx The index at which the frequency data is to be fetched.
     * @return The data matrix at the target index.
     */
    Eigen::MatrixXcd get_cplxData_at_f( unsigned int f_idx ) const;

    /**
     * Replace the block of continuous data within the frequency data array with the given block.
     * 
     * @param lead the starting index of the block.
     * @param newBlk_re insert block real part.
     * @param newBlk_im insert block imag part.
     */
    void set_cplxData_block( unsigned int lead, const Matrix3DXd& newBlk_re, const Matrix3DXd& newBlk_im );

    /**
     * Obtain the frequency vector.
     * 
     * @return The frequency vector (copy) of the current fData instance.
     */
    Eigen::VectorXd getF_vec() const;
    /**
     * Obtain the frequency data part A (Real part or magnitude).
     * 
     * @return The portion A (Real part or magnitude) of the frequency data (copy).
     */
    Matrix3DXd getXr_vec() const;
    /**
     * Obtain the frequency data part B (Imaginary part or phase).
     * 
     * @return The portion B (Imaginary part or phase) of the frequency data (copy).
     */
    Matrix3DXd getXi_vec() const;

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Printing
// ====================================================================== >>>>>

    /**
     * Write the data content of the current fData instance to the target file (.txt).
     * 
     * The format of the final data file follows a minimal version of a touchstone file.
     * For example, for a 3-by-3 data matrix case, the first column is the frequency column
     * and the subsequent columns follow the data order: 
     * 
     * X11a X11b X12a X12b X13a X13b X21a X21b etc.
     * 
     * @param fileDir The directory where the data file is created.
     * @param fileStem The file stem (file name without extension and directory) of the file to be created.
     * @param options Does nothing for now ...
     * 
     * @note 2-ports system follow the same data column ordering, and not the X11 X21 X22 X12 ordering.
     */
    void print_to( const string& fileDir, const string& fileStem, int options );

// ====================================================================== <<<<<

    /**
     * Serialize method to save current object state to a file.
     * 
     * @param fullFilename The full filename of the target file to create where 
     * the serialized data is to be stored.
     */
    void serialize(const std::string& fullFilename) const;

    /**
     * Serialize method using existing binary file stream.
     * 
     * @param ofs The binary file stream through which the serialized data of 
     * this fData instance is to be send.
     */
    void serialize( std::ofstream& ofs ) const;

    /**
     * Deserialize method to imprint the current fData instance state with existing 
     * state extracted from a binary file.
     * 
     * @param fullFilename The full filename of the target file from which the serialized
     * data is to be extracted from.
     */
    void deserialize( const std::string& fullFilename );

    /**
     * Deserialize method to imprint the current fData instance state with existing 
     * state extracted from an existing data file stream.
     * 
     * @param ifs The binary file stream through which a fData instance serialized data is to 
     * be extracted.
     */
    void deserialize( std::ifstream& ifs );

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

    // The vector of frequencies in hertz (unit saved separately).
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




