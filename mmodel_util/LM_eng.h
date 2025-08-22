#ifndef LM_ENG_H
#define LM_ENG_H

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>   
#include <vector> 

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

// ====================================================================== >>>>>
//      Static Support Functions
// ====================================================================== >>>>>

    /*
    Highly specific support function which creates a LM_eng based on the 
    data file specified by "fullFileName" and goes through the SFML process
    to generate singular values of real Loewner Matrix system.
    The singular values are then saved at the target directory "destDir" as
    a simple text file.
    */
    static shared_ptr<LM_eng> print_singVals( const string& fullFileName, 
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

    void step1_fData_partition();

    void step2_LM_construct();

    void step3_LM_re_trans();

    void step4_LM_pencil_SVD();

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
void set_fData( const fData& inData );

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
Eigen::MatrixXcd get_LM_re() const;
// Obtain the shifted-Loewner Matrix after real transform.
Eigen::MatrixXcd get_SLM_re() const;
// Obtain the partition 1 data row vector after real transform.
Eigen::MatrixXcd get_W_re() const;
// Obtain the partition 2 data column vector after real transform.
Eigen::MatrixXcd get_F_re() const;


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

    // The reduced frequency array size (training set size).
    unsigned int s1_fr_len = 100;
    // Frequency partition 1 has DC point.
    bool f1_has_DC_pt = false;
    // Frequency partition 2 has DC point.
    bool f2_has_DC_pt = false;
    // The reduced frequency set.
    shared_ptr<fData> myFr;
    // The two frequency data partitions.
    vector< unsigned int > partit1IdxArr;
    vector< unsigned int > partit2IdxArr;

};


namespace LM_UTIL{


    /*
    Construct a Loewner Matrix.
    */
    shared_ptr<Eigen::MatrixXcd> build_LM( const fData& f1Data, const fData& f2Data );

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