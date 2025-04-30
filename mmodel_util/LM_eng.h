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


    Eigen::MatrixXd build_LM( const fData& f1Data, const fData& f2Data );

}



#endif  // LM_ENG_H