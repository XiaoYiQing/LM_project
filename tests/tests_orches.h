#ifndef TESTS_ORCHESTRATORS_H
#define TESTS_ORCHESTRATORS_H




#include "orchestrators.h"

extern string RES_PATH_XYQ_str;
extern string SRC_PATH_XYQ_str;

using namespace std;


namespace tests{

    void orch_SFLM_full_run_test();

    void orch_SFLM_direc_re_run_test();

    void orch_singVal_extract_run_test();

}



#endif      // TESTS_ORCHESTRATORS_H