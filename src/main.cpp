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


/*
Function to retrieve frequency data from files ending with the extension of format
".sXp" where X is any positive integer representing number of I/O ports.
*/
void read_sXp_file( const string& fullFileName );


int main() {

    int err_ret_val = 1;

    cout << RES_PATH_XYQ_str << endl;

    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";


    read_sXp_file( fullFileName );


    return 0; // Indicates successful completion of the program

}



void read_sXp_file( const string& fullFileName ){

// ---------------------------------------------------------------------- >>>>>
//      Full Name Parsing
// ---------------------------------------------------------------------- >>>>>

    std::filesystem::path fullFilePath(fullFileName);
    std::string fileDir = fullFilePath.parent_path().string();
    std::string fileStem = fullFilePath.stem().string();
    std::string fileExt = fullFilePath.extension().string();
    
// ---------------------------------------------------------------------- <<<<<

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
       return;
   }

   
   int port_cnt_tmp = -1;
   try {
       port_cnt_tmp = std::stoi( xValue );
       const int loool = std::stoi( xValue );
       cout << "The integer is: " << port_cnt_tmp << std::endl;
   } catch (const std::invalid_argument& e) {
       cout << "Invalid input: The string does not contain a valid integer." << std::endl;
       cout << e.what() << endl;
       return;
   } catch (const std::out_of_range& e) {
       cout << "Invalid input: The integer is out of range." << std::endl;
       cout << e.what() << endl;
       return;
   }
   
   // Create the constant version of the port count.
   const unsigned int port_cnt = port_cnt_tmp;

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      File Parameters Read
// ---------------------------------------------------------------------- >>>>>

    // Open the input file stream.
    std::ifstream inputFile( fullFilePath );
    if( !inputFile ){
        cout << "Failed to read file." << endl;
        return;
    }else{
        cout << "File read successful." << endl;
    }


    // If the first character of a line is == comm_mark, the line is a comment.
    string comm_mark = "!";
    // If the first character of a line is == opt_line_mark, the line holds the various
    // options of the data format.
    string opt_line_mark = "#";
    // The variables holding the line and word currently read, respectively.
    string line, word;

    // Define the data options' default values.
    vector<string> options = { "GHZ", "S", "MA", "R", "50" };

    // Boolean flag for indicating whether the line stream has reached the data lines.
    bool data_reached = false;

    while( !data_reached && getline( inputFile, line ) ){

        // Set the stream for the current line.
        istringstream iss(line);
        // Read the first word.
        iss >> word;

        // Check for comment line mark.
        if( word == comm_mark ){
            cout << "This is a comment!" << endl;

        // Check for option line mark.
        }else if( word == opt_line_mark ){
            cout << "This is an option!" << endl;
            int opt_idx = 0;
            while( iss >> word ){
                options.at(opt_idx) = word;
                opt_idx++;
            }
            cout << endl;

        // Check for end of consecutive series of comment and option lines.
        }else{

            data_reached = true;

        }

    }

    // If we reached the end of the file without reaching any data line, abort.
    if( !data_reached ){
        cout << "The entire file has been read without reaching a data line." << endl;
        return;
    }

// ---------------------------------------------------------------------- <<<<<

    

// ---------------------------------------------------------------------- >>>>>
//      File Data Read
// ---------------------------------------------------------------------- >>>>>

    // Total number of parameter within the data matrix 
    // (For example, 2 by 2 S-parameters matrix has 4 individual S-parameters).
    unsigned int mat_ent_cnt = port_cnt*port_cnt;

    // File parsing control variables.
    unsigned int line_idx = 0;
    unsigned int data_idx = 0;
    unsigned int res_blk_size = 200;
    unsigned int curr_vec_size = 0;
    // Temporary data value to be used during translation from string to double.
    double tmp_val = 0;

    // The frequency vector.
    vector< double > f_vec;
    // The portion A data vector (A is typically magnitude or real part).
    vector< vector<double> > val_M_vec;
    // The portion B data vector (B is typically phase or imaginary part).
    vector< vector<double> > val_P_vec;

    f_vec.reserve( res_blk_size );
    val_M_vec.reserve( res_blk_size );
    val_P_vec.reserve( res_blk_size );

    // Initialize inner vectors and resize them to the desired size
    for (size_t i = 0; i < res_blk_size; i++) {
        val_M_vec.emplace_back( mat_ent_cnt );
        val_P_vec.emplace_back( mat_ent_cnt );
    }
    // Update vector size.
    curr_vec_size = val_M_vec.size();

    do{

        // Set the stream for the current line.
        istringstream iss( line );

        // Read the frequency word.
        iss >> word;
        // Translate the word into a double value freq.
        tmp_val = std::stod( word );
        // Save the frequency value.
        f_vec.push_back( tmp_val );


        for( unsigned int z = 0; z < mat_ent_cnt; z++ ){
            
            // Read the next data mag.
            iss >> word;    tmp_val = std::stod( word );
            val_M_vec.at( line_idx ).at( z ) = tmp_val;
            // Read the next data phase.
            iss >> word;    tmp_val = std::stod( word );
            val_P_vec.at( line_idx ).at( z ) = tmp_val;

        }

        line_idx++;
        if( line_idx >= curr_vec_size ){

            f_vec.reserve( line_idx + res_blk_size );
            val_M_vec.reserve( line_idx + res_blk_size );
            val_P_vec.reserve( line_idx + res_blk_size );

            for (size_t i = line_idx; i < line_idx + res_blk_size; i++) {
                val_M_vec.emplace_back( mat_ent_cnt );
                val_P_vec.emplace_back( mat_ent_cnt );
            }

            curr_vec_size += res_blk_size;

        }

    }while( getline( inputFile, line ) );
        
    // Deallocate unused reserved memory from the vectors.
    f_vec.shrink_to_fit();
    val_M_vec.resize( line_idx );
    val_M_vec.shrink_to_fit();
    val_P_vec.resize( line_idx );
    val_P_vec.shrink_to_fit();

// ---------------------------------------------------------------------- <<<<<

    return;

}