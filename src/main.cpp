// Preprocessor variable created when running CMAKE process. This variable 
// indicates the location where non-compilabe or non c/c++ support files are located.
#ifndef RES_PATH_XYQ
#   define RES_PATH_XYQ ""
#else
#endif

#include <iostream>
#include <filesystem>
#include <fstream>
#include <regex>
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

    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";

    std::filesystem::path fullFilePath(fullFileName);
    std::string fileDir = fullFilePath.parent_path().string();
    std::string fileStem = fullFilePath.stem().string();
    std::string fileExt = fullFilePath.extension().string();


// ---------------------------------------------------------------------- >>>>>
//      Port Count Regex Determination
// ---------------------------------------------------------------------- >>>>>
    /*
    Regex for determining the positive integer value X in the pattern ".sXp".
    */
    regex pattern(R"(\.s(\d+)p)");
    // The match result variable.
    smatch matches;
    // The exact number of the match in string.
    string xValue;

    if ( regex_match(fileExt, matches, pattern) ) {
        if (matches.size() > 1) { // Check if we have a match for the integer
            xValue = matches[1]; // Get the captured group (X)
            cout << "X: " << xValue << endl; // Output the value of X
        }
    } else {
        cout << "No match found." << endl;
        return 1;
    }

    int port_cnt = -1;
    try {
        port_cnt = std::stoi( xValue );
        std::cout << "The integer is: " << port_cnt << std::endl;
    } catch (const std::invalid_argument& e) {
        std::cout << "Invalid input: The string does not contain a valid integer." << std::endl;
    } catch (const std::out_of_range& e) {
        std::cout << "Invalid input: The integer is out of range." << std::endl;
    }
    
// ---------------------------------------------------------------------- <<<<<

    // Open the input file stream.
    std::ifstream inputFile( fullFilePath );
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

    // Define the data options' default values.
    vector<string> options = { "GHZ", "S", "MA", "R", "50" };

    // Boolean flag for indicating whether the line stream has reached the data lines.
    bool data_reached = false;

    while( !data_reached && std::getline( inputFile, line ) ){

        // Obtain the next word.
        std::string word;
        std::istringstream iss(line);

        iss >> word;
        // Check for comment line mark.
        if( word == comm_mark ){
            cout << "This is a comment!" << endl;
        }else if( word == opt_line_mark ){
            cout << "This is an option!" << endl;
            int opt_idx = 0;
            while( iss >> word ){
                options.at(opt_idx) = word;
                opt_idx++;
            }
            cout << endl;
        }else{
            data_reached = true;
        }

    }





    return 0; // Indicates successful completion of the program

}