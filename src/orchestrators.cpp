
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


shared_ptr<LM_eng> FCT_SCR::SFLM_full_run( const fData& src_data, 
    vector<unsigned int> f_r_idx_vec, vector<unsigned int> f1_idx_vec,
    vector<unsigned int> f2_idx_vec ){
    
    shared_ptr<LM_eng> my_LM_eng;
    
    // Initialization.
    try{
        my_LM_eng = std::make_shared<LM_eng>( LM_eng() );
    }catch( const std::exception& e ){
        cerr << "SFML process interrupted at engine creation: " << e.what() << endl;
    }

    // Step 0: data insertion.
    try{
        my_LM_eng->step0_fData_set( src_data, f_r_idx_vec );
    }catch( const std::invalid_argument& e ){
        cerr << "SFML process interrupted at step 0: " << e.what() << endl;
        return my_LM_eng;
    }catch( const std::out_of_range& e ){
        cerr << "SFML process interrupted at step 0: " << e.what() << endl;
        return my_LM_eng;
    }

    // Step 1: data partitioning.
    try{
        my_LM_eng->step1_fData_partition( f1_idx_vec, f2_idx_vec );
    }catch( const std::runtime_error& e ){
        cerr << "SFML process interrupted at step 1: " << e.what() << endl;
        return my_LM_eng;
    }catch( const std::invalid_argument& e ){
        cerr << "SFML process interrupted at step 1: " << e.what() << endl;
        return my_LM_eng;
    }catch( const std::out_of_range& e ){
        cerr << "SFML process interrupted at step 1: " << e.what() << endl;
        return my_LM_eng;
    }

    // Step 2: LM construction.
    try{ 
        my_LM_eng->step2_LM_construct();
    }catch( const std::runtime_error& e ){
        cerr << "SFML process interrupted at step 2: " << e.what() << endl;
        return my_LM_eng;
    }catch( const std::invalid_argument& e ){
        cerr << "SFML process interrupted at step 2: " << e.what() << endl;
        return my_LM_eng;
    }

    // Step 3: real matrix transform.
    try{
        my_LM_eng->step3_LM_re_trans();
    }catch( const std::runtime_error& e ){
        cerr << "SFML process interrupted at step 3: " << e.what() << endl;
        return my_LM_eng;
    }

    


    return my_LM_eng;

}
