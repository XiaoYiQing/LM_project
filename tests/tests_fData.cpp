#include "tests_fData.h"


using namespace std;

void tests::fData_test_1( unsigned int test_idx ){


    int case_cnt = 0;

    // 0- Base file read check.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // 0- Basic file reading test.
        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        cout << myF.get_f_cnt() << endl;
        cout << myF.get_f_scale_str() << endl;
        cout << myF.get_f_scale_num() << endl;
        cout << myF.get_reData_at_f( 10 ) << endl;
        cout << myF.get_imData_at_f( 10 ) << endl;

    }

    case_cnt++;
    // 1- Data Conversion check (DB <-> MA).
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // 0- Basic file reading test.
        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        cout << "[DB]" << endl;
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;

        // DB to MA
        cout << "[DB to MA]:" << endl;
        myF.data_format_Switch( fData::FDATA_FORMAT::MA );
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;

        // MA to DB
        cout << "[MA to DB]:" << endl;
        myF.data_format_Switch( fData::FDATA_FORMAT::DB );
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;
        

    }

    case_cnt++;
    // 2- Data Conversion check (DB <-> RI).
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // 0- Basic file reading test.
        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        cout << "[DB]" << endl;
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;

        // DB to RI
        cout << "[DB to RI]:" << endl;
        myF.data_format_Switch( fData::FDATA_FORMAT::RI );
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;

        cout << "[RI to DB]:" << endl;
        myF.data_format_Switch( fData::FDATA_FORMAT::DB );
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;

    }

    case_cnt++;
    // 3- Data Conversion check (MA <-> RI).
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // 0- Basic file reading test.
        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        cout << "[DB]" << endl;
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;

        // DB to RI
        cout << "[DB to RI]:" << endl;
        myF.data_format_Switch( fData::FDATA_FORMAT::RI );
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;

        cout << "[RI to MA]:" << endl;
        myF.data_format_Switch( fData::FDATA_FORMAT::MA );
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;

        cout << "[MA to RI]:" << endl;
        myF.data_format_Switch( fData::FDATA_FORMAT::RI );
        cout << myF.get_reData_at_f(100) << endl;
        cout << myF.get_imData_at_f(100) << endl;


    }


    case_cnt++;
    // 4- Alternative constructor test.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // 0- Basic file reading test.
        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        // Eigen::VectorXd& f_vec, Matrix3DXd& Xr_vec, Matrix3DXd& Xi_vec
        Eigen::VectorXd f_vec = myF.getF_vec();
        Matrix3DXd Xr_vec = myF.getXr_vec();
        Matrix3DXd Xi_vec = myF.getXi_vec();
        fData myF_copy = fData( f_vec, Xr_vec, Xi_vec );

        cout << myF.get_fval_at(51) << endl;
        cout << myF_copy.get_fval_at(51) << endl;

        bool match_bool = true;

        // Generate 10 random indices.
        vector<int>* myRandInt = randIntGen( 0, myF.get_f_cnt() - 1, 10 );

        // Compare the copy and the original data at these random sampling points.
        for( int z : *myRandInt ){
            match_bool = match_bool && ( myF.get_fval_at(z) == myF_copy.get_fval_at(z) );
            match_bool = match_bool && ( myF.get_reData_at_f(z) == myF_copy.get_reData_at_f(z) );
            match_bool = match_bool && ( myF.get_imData_at_f(z) == myF_copy.get_imData_at_f(z) );

            if( !match_bool ){
                break;
            }
        }
        cout << "Match result: " << match_bool << endl;

    }


    case_cnt++;
    // 5- get_cplxData_at_f() test.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        unsigned testIdx = 2;

        Eigen::MatrixXcd tarMat = myF.get_cplxData_at_f( testIdx );
        Eigen::MatrixXd tarMatRe = tarMat.real();
        Eigen::MatrixXd tarMatIm = tarMat.imag();

        bool match_bool = true;

        match_bool = match_bool && ( tarMatRe == myF.get_reData_at_f( testIdx ) );
        match_bool = match_bool && ( tarMatIm == myF.get_imData_at_f( testIdx ) );
        cout << "Complex data matching result: " << match_bool << endl;

        complex<double> tarF = myF.get_cplx_f_at( testIdx );
        cout << "Complex frequency: " << tarF << endl;
        fData myF_cconj = myF.gen_cplx_conj_set();
        cout << "Complex frequency: " << myF_cconj.get_cplx_f_at( testIdx ) << endl;

    }

    case_cnt++;
    // 6- data_prefix_switch test.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        cout << myF.get_fval_at( 10 ) << endl;
        myF.data_prefix_switch( fData::METRIC_PREFIX::G );
        cout << myF.get_fval_at( 10 ) << endl;


    }

}



void tests::fData_test_2( unsigned int test_idx ){

    int case_cnt = 0;

    // 0- Complex conjugate set generation
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );
        
        myF.gen_cplx_conj_set();

    }

    case_cnt++;
    // 1- Default partition generation
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        // Generate two partitions based on linear distribution.
        vector< shared_ptr<fData> > myPartits = myF.gen_2_partit();

        shared_ptr<fData> partit1 = myPartits.at(0);
        shared_ptr<fData> partit2 = myPartits.at(1);

        bool match_bool = true;
        vector< unsigned int > even_idx_arr = utils::gen_even_idx_arr( 0, myF.get_f_cnt() - 1 );
        vector< unsigned int > odd_idx_arr = utils::gen_odd_idx_arr( 0, myF.get_f_cnt() - 1 );

        for( unsigned int z = 0; z < even_idx_arr.size(); z++ ){
            match_bool = match_bool && ( partit1->get_fval_at(z) == myF.get_fval_at( even_idx_arr.at(z) ) );
            match_bool = match_bool && ( partit1->get_reData_at_f(z) == myF.get_reData_at_f( even_idx_arr.at(z) ) );
            match_bool = match_bool && ( partit1->get_imData_at_f(z) == myF.get_imData_at_f( even_idx_arr.at(z) ) );
        }
        for( unsigned int z = 0; z < odd_idx_arr.size(); z++ ){
            match_bool = match_bool && ( partit2->get_fval_at(z) == myF.get_fval_at( odd_idx_arr.at(z) ) );
            match_bool = match_bool && ( partit2->get_reData_at_f(z) == myF.get_reData_at_f( odd_idx_arr.at(z) ) );
            match_bool = match_bool && ( partit2->get_imData_at_f(z) == myF.get_imData_at_f( odd_idx_arr.at(z) ) );
        }
        
        cout << "Interleaving partitioned data match: " << match_bool << endl;

    }

    case_cnt++;
    // 2- Reduced freq. data set generation.
    if( test_idx == case_cnt ){

        // Subset size.
        unsigned int rSize = 100;

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        
        // Generate a linear index vector.
        vector< unsigned int > fr_idx_vec = utils::gen_lin_idx_arr( 0, myF.get_f_cnt() - 1, rSize );

        // Generate a reduced fData.
        shared_ptr<fData> myFr = myF.red_partit( fr_idx_vec );

        bool match_bool = true;
        for( unsigned int z = 0; z < fr_idx_vec.size(); z++ ){
            match_bool = match_bool && ( myFr->get_fval_at(z) == myF.get_fval_at( fr_idx_vec.at(z) ) );
            match_bool = match_bool && ( myFr->get_reData_at_f(z) == myF.get_reData_at_f( fr_idx_vec.at(z) ) );
            match_bool = match_bool && ( myFr->get_imData_at_f(z) == myF.get_imData_at_f( fr_idx_vec.at(z) ) );
        }
        cout << "Reduced frequency data object match: " << match_bool << endl;

    }


    case_cnt++;
    // 3- gen_cplx_conj_comb() test.
    if( test_idx == case_cnt ){

        // Subset size.
        unsigned int rSize = 100;

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );
        
        // Generate a linear index vector.
        vector< unsigned int > fr_idx_vec = utils::gen_lin_idx_arr( 0, myF.get_f_cnt() - 1, rSize );

        // Generate a reduced fData.
        shared_ptr<fData> myFr = myF.red_partit( fr_idx_vec );

        // Generate the complex conjugate set and insert into the reduced array.
        shared_ptr<fData> myFr_cconj = myFr->gen_cplx_conj_comb();

        // Define the test boolean result.
        bool match_bool = true;
        
        // Obtain the original reduced fData size.
        unsigned int orig_fr_cnt = myFr->get_f_cnt();

        unsigned int idx_offset = 0;
        if( myFr->hasDC() ){
            idx_offset = 1;
        }
        // Verify the placement of the data and their respective complex conjugates.
        for( unsigned int z = idx_offset; z < orig_fr_cnt; z++ ){

            if( !match_bool ){
                break;
            }

            unsigned int z2 = 2*z - idx_offset;

            match_bool = match_bool && ( myFr->get_fval_at(z) ==
                myFr_cconj->get_fval_at( z2 ) );
            match_bool = match_bool && ( myFr->get_fval_at(z) ==
                -1*myFr_cconj->get_fval_at( z2 + 1 ) );
            match_bool = match_bool && ( myFr->get_cplxData_at_f(z) ==  
                myFr_cconj->get_cplxData_at_f( z2 ) );
            match_bool = match_bool && ( myFr->get_cplxData_at_f(z) ==  
                ( myFr_cconj->get_cplxData_at_f( z2 + 1 ) ).conjugate() );

        }



        cout << "Self complex conjugate injected fData test: " << match_bool << endl;

    }

}


void tests::fData_enum_test( unsigned int test_idx ){

    int case_cnt = 0;
    /*
    Test get_METRIC_PREFIX_next.
    */
    if( test_idx == case_cnt ){

        // Initialize prefix variable.
        fData::METRIC_PREFIX nextPrefix = fData::METRIC_PREFIX::NONE;
        // Initialize test boolean.
        bool test_bool = true;

        // Standard test.
        nextPrefix = fData::get_METRIC_PREFIX_next( fData::METRIC_PREFIX::k, true );
        test_bool = test_bool && nextPrefix == fData::METRIC_PREFIX::M;
        nextPrefix = fData::get_METRIC_PREFIX_next( fData::METRIC_PREFIX::k, false );
        test_bool = test_bool && nextPrefix == fData::METRIC_PREFIX::h;
        if( test_bool ){
            cout << "get_METRIC_PREFIX_next standard test: passed!" << endl;
        }else{
            cout << "get_METRIC_PREFIX_next standard test: failed!" << endl;
        }

        // Reset test boolean.
        test_bool = true;
        // Try to get next higher prefix at the highest defined prefix enum.
        fData::METRIC_PREFIX testPrefix = static_cast<fData::METRIC_PREFIX>( fData::METRIC_PREFIX_Count - 1 );
        try{
            nextPrefix = fData::get_METRIC_PREFIX_next( testPrefix, true );
        }catch( const std::out_of_range& e ){
            cout << e.what() << endl;
            test_bool = test_bool && true;
        }catch( ... ){
            test_bool = test_bool && false;
        }
        
        // Try to get next lower prefix at the highest defined prefix enum.
        testPrefix = static_cast<fData::METRIC_PREFIX>( 0 );
        try{
            nextPrefix = fData::get_METRIC_PREFIX_next( testPrefix, false );
        }catch( const std::out_of_range& e ){
            cout << e.what() << endl;
            test_bool = test_bool && true;
        }catch( ... ){
            test_bool = test_bool && false;
        }

        if( test_bool ){
            cout << "get_METRIC_PREFIX_next out of bound failure test: passed!" << endl;
        }else{
            cout << "get_METRIC_PREFIX_next out of bound failure test: failed!" << endl;
        }
        

    }

}


void tests::fData_test_sXp_read( unsigned int test_idx ){

    int case_cnt = 0;

    // 2-port test case.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/test_res_dir/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        bool test1_bool = true;

        test1_bool = test1_bool && myF.get_f_cnt() == 500;
        test1_bool = test1_bool && myF.get_f_scale_str() == 
            fData::get_METRIC_PREFIX_Str( fData::METRIC_PREFIX::NONE );
        test1_bool = test1_bool && myF.get_f_scale_num() == 1;

        Eigen::MatrixXd mag1(2,2);
        mag1 << -14.56544180548246, -0.1650035838894745,
            -0.1650035838894637, -14.56394882982146;
        test1_bool = test1_bool && myF.get_reData_at_f( 10 ) == mag1;

        Eigen::MatrixXd phase1(2,2);
        phase1 << 18.51450720004394, -70.84049885284401,
                    -70.84049885284401, 20.32410950236233;
        double myPI = 2*std::asin(1.0);
        phase1 = phase1 *( myPI/180 );
        test1_bool = test1_bool && 
            ( myF.get_imData_at_f( 10 ) - phase1 ).cwiseAbs().maxCoeff() < 1e-9;

        if( test1_bool ){
            cout << "Test 1 [.s2p file reading] passed!" << endl;
        }else{
            cout << "Test 1 [.s2p file reading] failed!" << endl;
        }
    
    }

    case_cnt++;
    // 4-port test case.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        string fullFileName = RES_PATH_XYQ_str + "/bondwire_with_strip_design3.s4p";
        fData::read_sXp_file( myF, fullFileName );

        bool test_bool = true;

        test_bool = test_bool && myF.get_f_cnt() == 500;
        test_bool = test_bool && myF.get_f_scale_str() == 
            fData::get_METRIC_PREFIX_Str( fData::METRIC_PREFIX::G );
        test_bool = test_bool && myF.get_f_scale_num() == 1e9;

        Eigen::MatrixXd mag1(4,4);
        mag1 << 0.958817972108667, 0.127454551904343, 0.0019875694004973, 0.00198918738171939, 
            0.127454551904343, 0.958795780807044, 0.0019278729865623, 0.00192596500436221,
            0.00198756940049718, 0.00192787298656213, 0.0390404547413966, 0.983742772847083, 
            0.00198918738171929, 0.00192596500436209, 0.983742772847081, 0.039015297765394;
        test_bool = test_bool && myF.get_reData_at_f( 10 ) == mag1;

        Eigen::MatrixXd phase1(4,4);
        phase1 << -51.8369933305439, 29.0982937219069, 42.4317772153111, 42.5781673192973, 
                29.098293721907,  -52.1954463168092, 42.6420711488167, 42.4834166003759, 
                42.4317772153035, 42.6420711488151, -147.022259604695, -54.324770082114, 
                42.5781673192937, 42.4834166003802, -54.3247700821141, -147.059096517886;
        double myPI = 2*std::asin(1.0);
        phase1 = phase1 *( myPI/180 );
        test_bool = test_bool && 
            ( myF.get_imData_at_f( 10 ) - phase1 ).cwiseAbs().maxCoeff() < 1e-9;

        if( test_bool ){
            cout << "Test 2 [.s4p file reading] passed!" << endl;
        }else{
            cout << "Test 2 [.s4p file reading] failed!" << endl;
        }

    }

}


void tests::fData_setFunc_tests( unsigned int test_idx ){

    int case_cnt = 0;

    // 0- Test set IO.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        unsigned int new_out = 4;
        unsigned int new_in = 3;
        myF.set_IO_cnt( new_in, new_out );

        bool match_bool = true;
        match_bool = match_bool && ( myF.get_f_cnt() == 0 );
        match_bool = match_bool && ( myF.get_out_cnt() == new_out );
        match_bool = match_bool && ( myF.get_in_cnt() == new_in );
        cout << "Set IO expected outcome match: " << match_bool << endl;

    }

    case_cnt++;
    // 1- Test set_fval_block.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        // Set the block's insert start point.
        unsigned int lead = 10;
        // Set the block size (Number of frequency data to insert).
        unsigned int block_size = 4;

        Eigen::VectorXd f_blk = Eigen::VectorXd( block_size );
        for( unsigned int z = 0; z < block_size; z++ ){
            f_blk(z) = z * utils::rDoubleGen( 0.0, 10.0, 1 )->at(0);
        }

        myF.set_fval_block( lead, f_blk );

        bool match_bool = true;
        for( unsigned int z = 0; z < block_size; z++ ){
            
            match_bool = match_bool && ( myF.get_fval_at(z) == f_blk(z) );
            if( !match_bool ){
                break;
            }

        }

        cout << "Freq Block insert test match: " << match_bool << endl;

    }


    case_cnt++;
    // 2- Test set_cplxData_block.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        unsigned int row_cnt = myF.get_out_cnt();
        unsigned int col_cnt = myF.get_in_cnt();

        

        // Set the block's insert start point.
        unsigned int lead = 10;
        // Set the block size (Number of frequency data to insert).
        unsigned int block_size = 4;

        Matrix3DXd newBlk_re = Matrix3DXd( row_cnt, col_cnt, block_size );
        Matrix3DXd newBlk_im = Matrix3DXd( row_cnt, col_cnt, block_size );
        for( unsigned int z = 0; z < block_size; z++ ){
            Eigen::MatrixXd mat_z = Eigen::MatrixXd( row_cnt, col_cnt );
            mat_z.setOnes();
            mat_z *= z * utils::rDoubleGen( 0.0, 10.0, 1 )->at(0);
            newBlk_re.set( z, mat_z );
            newBlk_im.set( z, mat_z );
        }

        myF.set_cplxData_block( lead, newBlk_re, newBlk_im );


        bool match_bool = true;
        for( unsigned int z = 0; z < block_size; z++ ){
            
            Eigen::MatrixXcd mat_z = myF.get_cplxData_at_f( z );
            match_bool = match_bool && ( mat_z.real() == newBlk_re.at(z) );
            match_bool = match_bool && ( mat_z.imag() == newBlk_im.at(z) );
            if( !match_bool ){
                break;
            }

        }

        cout << "Freq Data Block insert test match: " << match_bool << endl;

    }
    

    case_cnt++;

}




void tests::fData_LTspice_data_read_test(){

    // Define our frequency data object.
    fData myF;

    // 0- Basic file reading test.
    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/audioamp/audioamp.txt";
    fData::read_LTspice_Sp_file( myF, fullFileName );

    cout << myF.get_f_cnt() << endl;
    cout << myF.get_f_scale_str() << endl;
    cout << myF.get_f_scale_num() << endl;

    cout << "dB data:" << endl;
    cout << myF.get_reData_at_f( 10 ) << endl;
    cout << myF.get_imData_at_f( 10 ) << endl;

    cout << "RI data:" << endl;
    myF.data_format_Switch( fData::FDATA_FORMAT::RI );  
    cout << myF.get_reData_at_f( 10 ) << endl;
    cout << myF.get_imData_at_f( 10 ) << endl;
    

}


void tests::fData_print_test( unsigned int test_idx ){

    int case_cnt = 0;
    // 0- Standard touchstone test case.
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;

        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
        fData::read_sXp_file( myF, fullFileName );

        string targetDir = "C:/Users/Yi Qing Xiao/Documents/Cpp_projects/LM_project/data_output";
        string targetStemName = "tmp_data_file";

        myF.print_to( targetDir, targetStemName, 0 );

    }

}   


void tests::fData_prefix_manip_test( unsigned int test_idx ){

    int case_cnt = 0;
    // 0- get_METRIC_PREFIX_for_val tests
    if( test_idx == case_cnt ){

        bool test_bool = true;
        fData::METRIC_PREFIX curr_pref;
        string err_msg = "get_METRIC_PREFIX_for_val test: failed at ";

        char buf[64];

        vector<double> testVal_vec = { 1e-16, 1e-12, 1e-9, 1e-6, 
            1e-3, 1e-2, 1e-1, 1, 1e+1, 1e2, 1e3, 1e6, 1e9, 1e12 };
        vector<fData::METRIC_PREFIX> test_pref_vec = {
            fData::METRIC_PREFIX::p,
            fData::METRIC_PREFIX::p,
            fData::METRIC_PREFIX::n,
            fData::METRIC_PREFIX::mu,
            fData::METRIC_PREFIX::m,
            fData::METRIC_PREFIX::c,
            fData::METRIC_PREFIX::d,
            fData::METRIC_PREFIX::NONE,
            fData::METRIC_PREFIX::da,
            fData::METRIC_PREFIX::h,
            fData::METRIC_PREFIX::k,
            fData::METRIC_PREFIX::M,
            fData::METRIC_PREFIX::G,
            fData::METRIC_PREFIX::T
        };

        for( unsigned int z = 0; z < testVal_vec.size(); z++ ){
            curr_pref = fData::get_METRIC_PREFIX_for_val( testVal_vec.at(z) );
            if( curr_pref != test_pref_vec.at(z) ){
                test_bool = false;
                std::snprintf( buf, sizeof(buf), "%.6e", testVal_vec.at(z) );
                cout << err_msg + string( buf ) << endl;
            }
        }

        if( test_bool ){
            cout << "get_METRIC_PREFIX_for_val test: passed!" << endl;
        }

    }


    // 1- f_normalize tests
    case_cnt++;
    if( test_idx == case_cnt ){

        // Define our frequency data object.
        fData myF;
        // Define the full file name.
        string fullFileName = RES_PATH_XYQ_str + "/audioamp/audioamp.txt";
        // Parse target data file and insert freq data into the fData object.
        fData::read_LTspice_Sp_file( myF, fullFileName );
        
        // Switch the data format to real/imaginary.
        myF.data_format_Switch( fData::FDATA_FORMAT::RI );

        // Perform normalization of the freq vector.
        myF.f_normalize();

    }

    


    

}