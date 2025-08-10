#ifndef TEST_NUMUTILS_H
#define TEST_NUMUTILS_H


#include <iostream>
#include <string>

#include "numUtils.h"
#include "eigenUtils.h"

extern string RES_PATH_XYQ_str;

using namespace std;


namespace tests{

    void numUtils_test_1( unsigned int test_idx );

    // Test functions in the eigenUtils file.
    void eigenUtils_test_1( unsigned int test_idx );

}





#endif  // TEST_NUMUTILS_H