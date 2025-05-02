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
    Construct the row matrix vector containing the partition 2 data matrices in 
    the same order they are used to construct the LM and SLM.
    */
    shared_ptr<Eigen::MatrixXcd> build_W( const fData& f2Data );

}



#endif  // LM_ENG_H