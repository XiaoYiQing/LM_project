#ifndef LM_ENG_H
#define LM_ENG_H

#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>   
#include <vector> 

#include "fData.h"


using namespace std;


class LM_eng{

public: 

    LM_eng();

protected:

    // The starting frequency data from which the LM is to be constructed.
    fData orig_fData;

    

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
    

}




#endif  // LM_ENG_H