#ifndef ORCHESTRATORS_H
#define ORCHESTRATORS_H



#include <string>
#include <vector>

#include "fData.h"
#include "LM_eng.h"


using namespace std;

extern string RES_PATH_XYQ_str;
extern string SRC_PATH_XYQ_str;


namespace FCT_SCR{

    /*
    Function follows a predefined list of S-param data files and creates a
    LM_eng for each file and then print the real LM pencil singular values to
    a predefined data output directory.
    */
    void singVal_extract_run();



};



#endif  // ORCHESTRATORS_H