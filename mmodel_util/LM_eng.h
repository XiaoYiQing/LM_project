#ifndef LM_ENG_H
#define LM_ENG_H

#include <algorithm>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>   
#include <vector> 

#include "eigenUtils.h"
#include "fData.h"
#include "LTI_descSyst.h"

using namespace std;


class LM_eng{

public: 

    /*
    Numerical threshold utilized by instances of this class for determining nullity.
    If a number's magnitude is smaller than this threshold, it is considered a zero.
    */
    static const double NUM_THRESH;

    /*
    Standard reduced frequency set size. This is the number of frequency entries to be
    extracted from the original fData set to create a subset fData. This subset fData
    is the one used to construct the SFLM system, because the full set fData is usually
    too large and overkill.
    */
    static const unsigned int STD_RED_FSET_SIZE;

// ====================================================================== >>>>>
//      Data Printing Functions
// ====================================================================== >>>>>

    /*
    Highly specific support function which creates a LM_eng based on the 
    data file specified by "fullFileName" and goes through the SFML process
    to generate singular values of real Loewner Matrix system.
    The singular values are then saved at the target directory "destDir" as
    a simple text file using the same file name stem as "fullFileName" with added
    suffix "_sv".
    For example, if 
        - destDir = "C:/project/singVal"
        - fullFileName = "lowPassFilterV2.s2p" or "lowPassFilterV2.txt"

    Then the generated singular values file will have the full file name:
        - "C:/project/singVal/lowPassFilterV2_vs.txt"
    */
    static shared_ptr<LM_eng> print_singVals( const string& fullFileName, 
        const string& destDir );
    
    /*
    Print the singular values of the target LM system at the specific directory under
    the filename stem.
    For example, if 
        - destDir = "C:/project/singVal"
        - fileStem = "lowPassFilterV2"

    Then the generated singular values file will have the full file name:
        - "C:/project/singVal/lowPassFilterV2.txt"
    NOTE: The given "tar_LM_eng" must have generated the singular values already, or
    the process is aborted.
    */
    static void print_singVals( const LM_eng& tar_LM_eng, const string& fileStem, 
        const string& destDir );

// ====================================================================== <<<<<

// ====================================================================== >>>>>
//      Constructor
// ====================================================================== >>>>>

    LM_eng();

    LM_eng( const fData& inData );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Major LM System Steps
// ====================================================================== >>>>>

    /*
    The pre-SFML process step, which is simply setting the frequency data as well
    as selecting the portion of the data being used to construct the LMs and the rest
    for validation.
    */ 
    void step0_fData_set( const fData& inData );
    /*
    The pre-SFML process step, which is simply setting the frequency data as well
    as selecting the portion of the data being used to construct the LMs and the rest
    for validation.
    This function version allows user to specify which portion of the data is used
    to construct the LMs via specifying "fr_idx_arr_in".
    */ 
    void step0_fData_set( const fData& inData, const vector<unsigned int>& fr_idx_arr_in );

    /*
    SFML step where the data selected to construct the LMs are partitioning into 2 
    partitions.
    */
    void step1_fData_partition();
    /*
    SFML step where the data selected to construct the LMs are partitioning into 2 
    partitions.
    This function version lets the user directly dictate how the two partitions
    are created. Note that the two index arrays argument must:
        - have mutually exclusive indices.
        - have no repeated indices.
        - have a total of indices equal to the amount of data selected to construct the LMs.
    */
    void step1_fData_partition( const vector<unsigned int>& f1IdxVec, 
        const vector<unsigned int>& f2IdxVec );
    
    /*
    SFML step where the LM are constructed using the complex frequency data.
    As such, the LMs generated at this step are complex and not purely real.
    */
    void step2_LM_construct();
    
    /*
    SFLM step where the complex LMs are transformed into real LMs.
    */
    void step3_LM_re_trans();
    
    /*
    Special SFLM step where the real Loewner Matrices are constructed directly.
    However, the complex Loewner Matrices are skipped entirely, and are thus 
    inaccessible.
    */
    void step3skip2_LM_re_construct();
    
    /*
    SFLM step where the LM pencil is created and then undergoes singular value 
    decomposition.
    */
    void step4_LM_pencil_SVD();
    /*
    SFLM step where the LM pencil is created and then undergoes singular value 
    decomposition.
    This function version let's the user select the reference frequency.
    NOTE: the value would normally be chosen from the frequency array with which
    the LMs were constructed with.
    */
    void step4_LM_pencil_SVD( double f_ref );

    /*
    Given the number of singular values to be kept, create the transfer function
    of equal order from the LM system.
    NOTE: the input must be the number of singular values, NOT the singular value index
    where the cut-off occurs.
    */
    shared_ptr<LTI_descSyst> step5_LM_to_tf( unsigned int svd_ret_cnt );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Access Function
// ====================================================================== >>>>>

/*
Get the flag of the specified LM step index.
*/
bool get_flag( unsigned int flagIdx ) const;

// Insert new frequency data to the current LM engine.
// void set_fData( const fData& inData );

// Obtain the current fData.
fData get_fData() const;

/*
Obtain the reference frequency magnitude used to construct the LM pencil.
*/
double get_ref_f_mag() const;

// Obtain the number of outputs.
unsigned int get_out_cnt() const;
// Obtain the number of inputs.
unsigned int get_in_cnt() const;

// Obtain the frequency partition 1 data.
shared_ptr<fData> get_Fr1() const;
// Obtain the frequency partition 2 data.
shared_ptr<fData> get_Fr2() const;
// Obtain the complex conjugate frequency partition 1 data.
shared_ptr<fData> get_Frc1() const;
// Obtain the complex conjugate frequency partition 2 data.
shared_ptr<fData> get_Frc2() const;

// Obtain the Loewner Matrix.
Eigen::MatrixXcd get_LM() const;
// Obtain the shifted-Loewner Matrix.
Eigen::MatrixXcd get_SLM() const;
// Obtain the partition 1 data row vector.
Eigen::MatrixXcd get_W() const;
// Obtain the partition 2 data column vector.
Eigen::MatrixXcd get_F() const;

// Obtain the Loewner Matrix after real transform.
Eigen::MatrixXd get_LM_re() const;
// Obtain the shifted-Loewner Matrix after real transform.
Eigen::MatrixXd get_SLM_re() const;
// Obtain the partition 1 data row vector after real transform.
Eigen::MatrixXd get_W_re() const;
// Obtain the partition 2 data column vector after real transform.
Eigen::MatrixXd get_F_re() const;


// Obtain the current Loewner Matrix pencil's computed singular values.
Eigen::VectorXd get_singVals() const;
// Obtain the left singular vectors generated from the current LM pencil.
Eigen::MatrixXd get_U() const;
// Obtain the left singular vectors generated from the current LM pencil.
Eigen::MatrixXd get_V() const;

// ====================================================================== <<<<<


protected:

    bool flag0_data_set = false;
    bool flag1_data_prep = false;
    bool flag2_LM_const = false;
    bool flag3_re_trans = false;
    bool flag4_pen_SVD = false;

    // The Loewner Matrices.
    Eigen::MatrixXcd LM;
    Eigen::MatrixXcd SLM;
    Eigen::MatrixXcd W;
    Eigen::MatrixXcd F;
    // The Real Loewner Matrices.
    Eigen::MatrixXd LM_re;
    Eigen::MatrixXd SLM_re;
    Eigen::MatrixXd W_re;
    Eigen::MatrixXd F_re;
    // The reference frequency magnitude.
    double ref_f_mag = 0;
    // The singular values of the LM pencil.
    Eigen::VectorXd singVals;
    // The singular vectors of the LM pencil.
    Eigen::MatrixXd U, V;

    // The starting frequency data from which the LM is to be constructed.
    fData myFData;

    // Frequency partition 1 has DC point.
    bool f1_has_DC_pt = false;
    // Frequency partition 2 has DC point.
    bool f2_has_DC_pt = false;

    // The reduced frequency set.
    // shared_ptr<fData> myFr;
    // The reduced frequency set index array.
    vector< unsigned int > fr_idx_arr;
    /*
    The index array of partition 1 of the reduced frequency set.
    NOTE: Indexing with respect to the full fData myFData.
    */
    vector< unsigned int > partit1IdxArr;
    /*
    The index array of partition 2 of the reduced frequency set.
    NOTE: Indexing with respect to the full fData myFData.
    */
    vector< unsigned int > partit2IdxArr;

};


namespace LM_UTIL{

    /*
    Create two index arrays having mutually exclusive indices which
    can serve to create two partitions from the existing original set.
    */
    vector< vector< unsigned int > > gen_2_partit_idx_arr( unsigned int origSize );

    /*
    Construct a Loewner Matrix.
    */
    shared_ptr<Eigen::MatrixXcd> build_LM( const fData& f1Data, const fData& f2Data );

    /*
    Construct a real Loewner Matrix. A number of points need to be considered when
    utilizing this function:
    - The given frequency data is going to be inflated with the complete complex 
        conjugate set (except for the DC point, if present).
    - If the DC point is present, it must be the first entry of either f1Data or f2Data.
    - The real Loenwer Matrix is constructed directly rather than through transform 
        matrix multiplication.
    */
    shared_ptr<Eigen::MatrixXd> build_LM_re( const fData& f1Data, const fData& f2Data );

    /*
    Construct a shifted-Loewner Matrix.
    */
    shared_ptr<Eigen::MatrixXcd> build_SLM( const fData& f1Data, const fData& f2Data );

    /*
    Construct the row matrix vector containing the partition 1 data matrices in 
    the same order they are used to construct the LM and SLM.
    */
    shared_ptr<Eigen::MatrixXcd> build_W( const fData& f1Data );

    /*
    Construct the col matrix vector containing the partition 2 data matrices in 
    the same order they are used to construct the LM and SLM.
    */
    shared_ptr<Eigen::MatrixXcd> build_F( const fData& f2Data );

    /*
    Generate the Loenwer Matrix pencil.
    The pencil is constructed using the Loewner Matrix (LM), the shifted 
    Loewner Matrix (SLM), and a select reference frequency value that can be
    picked as any complex frequency that participated in contructing the LM.
    */
    shared_ptr<Eigen::MatrixXcd> build_LM_pencil( complex<double>, const Eigen::MatrixXcd& LM, 
        const Eigen::MatrixXcd& SLM );

    /*
    Generate the Loenwer Matrix pencil, but after the real transform.
    The pencil is constructed using the Loewner Matrix (LM), the shifted 
    Loewner Matrix (SLM), and a select reference frequency value that can be
    picked as any frequency magnitude that participated in contructing the LM.
    */
    shared_ptr<Eigen::MatrixXd> build_LM_pencil( double, const Eigen::MatrixXd& LM, 
        const Eigen::MatrixXd& SLM );

    /*
    Construct a real transformation matrix in the context of Loewner Matrix real transform.
    The transformation is based on a target vector of matrices, where sub_mat_size is defined
    as the length each matrix takes from the vector total length.
    Note that the transformation matrix is scaled by 2^(-0.5).
    The vector of matrices is expected to have the specific pair pattern where each 
    distinct matrix entry is immediately followed by its complex conjugate.
    In this case, the vector total length is found to be:
        total_len = 2*sub_cplx_blk_cnt*sub_mat_size.

    The only exception to the above pattern occurs when has_DC_pt = true, in which
    case the very first matrix entry of the vector is purely real. Naturally, it won't
    be followed by a complex conjugate version, so the vector total length would be:
        total_len = ( 2*sub_cplx_blk_cnt + 1 )*sub_mat_size.

    > has_DC_pt: If the DC point is present, it means the first sub block vector is real and need
        no complex to real transform.
    > sub_mat_size: The size of each sub-matrix.
    > sub_cplx_blk_cnt: The number of distinct sub-matrices, excluding the purely real DC point one, if present.
    */
    shared_ptr<Eigen::MatrixXcd> build_reT_mat( bool has_DC_pt, unsigned int sub_mat_size, unsigned int sub_cplx_blk_cnt );
    

};




#endif  // LM_ENG_H