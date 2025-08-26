
#include "orchestrators.h"





void FCT_SCR::singVal_extract_run(){

    // Define the destination directory.
    string data_dest_dir = SRC_PATH_XYQ_str + "/data_output/singVals_depo";

    // Create the stem, path, and ext arrays for the file.
    vector<string> file_stem_arr;
    vector<string> file_path_arr;
    vector<string> file_ext_arr;

    unsigned int file_cnt = 0;

    file_stem_arr.push_back( "Slink_a=100um_b=400um" );
    file_path_arr.push_back( RES_PATH_XYQ_str + "/" );
    file_ext_arr.push_back( ".s2p" );
    file_cnt++;

    file_stem_arr.push_back( "Slink_a=100um_b=425um" );
    file_path_arr.push_back( RES_PATH_XYQ_str + "/" );
    file_ext_arr.push_back( ".s2p" );
    file_cnt++;

    string fullFileName_z = "";
    for( unsigned int z = 0; z < file_cnt; z++ ){

        // Assemble current full filename.
        fullFileName_z = file_path_arr.at(z) + file_stem_arr.at(z) + file_ext_arr.at(z);

        // Perform the automated <data parsing> + <SFML process> + 
        // <singVals print> process.
        LM_eng::print_singVals( fullFileName_z, data_dest_dir );

    }


}
