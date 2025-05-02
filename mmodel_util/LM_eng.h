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
    Construct a real transformation matrix.
    > has_DC_pt: If the DC point is present, it means the first sub block vector is real and need
        no complex to real transform.
    > sub_blk_size: The number of rows (sub_blk_size.at(0)) and columns (sub_blk_size.at(1)) of a
        sub-matrix block.
    > sub_blk_cnt: The number of sub-matrix blocks.
    */
    shared_ptr<Eigen::MatrixXcd> build_reT_mat( bool has_DC_pt, vector<unsigned int> sub_blk_size, unsigned int sub_blk_cnt );
    

}



#endif  // LM_ENG_H