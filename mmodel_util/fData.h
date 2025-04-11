#ifndef FDATA_H
#define FDATA_H


#include <Eigen/Dense>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>   
#include <vector>  

using namespace std;

class fData{    

public:

    /*
    Function to retrieve frequency data from files ending with the extension of format
    ".sXp" where X is any positive integer representing number of I/O ports.
    */
    static void read_sXp_file( const string& fullFileName );

// ====================================================================== >>>>>
//      Constructors
// ====================================================================== >>>>>

    fData();

// ====================================================================== <<<<<

protected:



private:


};






#endif  // FDATA_H




