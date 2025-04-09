// Preprocessor variable created when running CMAKE process. This variable 
// indicates the location where non-compilabe or non c/c++ support files are located.
#ifndef RES_PATH_XYQ
#   define RES_PATH_XYQ ""
#else
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>   
#include <vector>  


#include "numUtils.h"


using namespace std;

// Global string variable holding the full path of the extra resource directory.
string RES_PATH_XYQ_str = string( RES_PATH_XYQ );


int main() {

    
    cout << RES_PATH_XYQ_str << endl;

    vector<int> tmp = randIntVectGen( 0, 5, 7 );

    string filename = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
    cout << "File name: " << filename << endl;

    std::ifstream inputFile( filename );
    if( !inputFile ){
        cout << "Failed to read file." << endl;
    }else{
        cout << "File read successful." << endl;
    }

    return 0; // Indicates successful completion of the program

}