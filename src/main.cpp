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


    // If the first character of a line is == comm_mark, the line is a comment.
    string comm_mark = "!";
    // If the first character of a line is == opt_line_mark, the line holds the various
    // options of the data format.
    string opt_line_mark = "#";

    // The variable holding the line currently read.
    string line;

    // Boolean flag for indicating whether the line stream has reached the data lines.
    bool data_reached = false;

    while( !data_reached && std::getline( inputFile, line ) ){

        // Obtain the next word.
        std::string word;
        std::istringstream iss(line);

        iss >> word;
        // Check for comment line mark.
        if( word == comm_mark ){
            cout << "\tThis is a comment!" << endl;
        }
        // Use getline to read words, separated by spaces
        while ( iss >> word ) {
            
        }

    }


    return 0; // Indicates successful completion of the program

}