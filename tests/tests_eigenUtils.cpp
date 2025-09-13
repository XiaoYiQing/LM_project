#include "tests_eigenUtils.h"


void tests::eigenUtils_test_1( unsigned int test_idx ){

    // Initialize test case index.
    int case_cnt = 0;

    int case_idx = 0;
    // 0- Test gen_lin_idx_arr.
    if( case_cnt == case_idx ){

        unsigned int num_cnt = 10;

        shared_ptr< vector<double> > tmp_vec = utils::rDoubleGen( -10.0, 10.0, num_cnt );

        // Map the std::vector data to Eigen::VectorXd
        Eigen::Map<Eigen::VectorXd> eigenVec(tmp_vec->data(), tmp_vec->size());

        // If you need an independent copy
        Eigen::VectorXd tarVec = eigenVec;

        string targetDir = SRC_PATH_XYQ_str + "/data_output";
        string targetStemName = "tmp_data_file";

        utils::VectorXd_to_file( targetDir, targetStemName, tarVec, 10 );

    }

}


void tests::gen_rand_MatrixXd_test(){

    Eigen::MatrixXd tmp = utils::gen_rand_MatrixXd( 4, 4 );
    cout << tmp << endl;

}


void tests::Vd_to_file_test(){

    unsigned int num_cnt = 10;

    vector<double> tmp_vec = *utils::rDoubleGen( -10.0, 10.0, num_cnt );

    string targetDir = SRC_PATH_XYQ_str + "/data_output";
    string targetStemName = "Vd_to_file_res";

    utils::Vd_to_file( targetDir, targetStemName, tmp_vec, 10 );
    
}


void tests::file_to_vec_test(){

// ---------------------------------------------------------------------- >>>>>
//      Initialization
// ---------------------------------------------------------------------- >>>>>

    string targetDir = RES_PATH_XYQ_str + "/test_res_dir";
    string targetStemName = "test_vec_data";
    string targetExt = ".txt";
    string fullFileName = targetDir + "/" + targetStemName + targetExt;

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Simple Correct Text File Read
// ---------------------------------------------------------------------- >>>>>

    

    vector<double> my_vec = utils::file_to_Vd( fullFileName );
    vector<double> my_vec_ans = {
        -1.4273852567e-01, -2.1484261800e+00, +9.6893777105e+00, -7.8754513961e-01,
        -8.3398261166e+00, +1.2544720764e+00, -3.8343728561e+00, -9.3529739505e-01,
        -9.7646242005e+00, -2.9892813255e+00
    };

    bool test_vect = my_vec.size() == my_vec_ans.size();
    for( unsigned int z = 0; z < my_vec.size(); z++ ){
        test_vect = test_vect && ( my_vec[z] == my_vec_ans[z] );
    }
    if( test_vect ){
        cout << "file_to_Vd correct read test: passed!" << endl;
    }else{
        cout << "file_to_Vd correct read test: failed!" << endl;
    }
    cout << endl;

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Simple Correct CSV File Read
// ---------------------------------------------------------------------- >>>>>

    targetExt = ".csv";
    fullFileName = targetDir + "/" + targetStemName + targetExt;

    vector<double> my_vec2 = utils::file_to_Vd( fullFileName );

    test_vect = my_vec2.size() == my_vec_ans.size();
    for( unsigned int z = 0; z < my_vec2.size(); z++ ){
        test_vect = test_vect && ( abs( my_vec2[z] - my_vec_ans[z] ) < 1e-12 );
    }
    if( test_vect ){
        cout << "file_to_Vd correct csv file read test: passed!" << endl;
    }else{
        cout << "file_to_Vd correct csv file read test: failed!" << endl;
    }
    cout << endl;

// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      Erroneous Text File Read
// ---------------------------------------------------------------------- >>>>>

    targetDir = RES_PATH_XYQ_str + "/test_res_dir";
    targetStemName = "test_bad_vec_data";
    targetExt = ".txt";

    fullFileName = targetDir + "/" + targetStemName + targetExt;
    try{
        vector<double> my_vec3 = utils::file_to_Vd( fullFileName );
        test_vect = false;
    }catch(...){
        test_vect = true;
    }
    if( test_vect ){
        cout << "file_to_Vd bad file read test: passed!" << endl;
    }else{
        cout << "file_to_Vd bad file read test: failed!" << endl;
    }
    cout << endl;

// ---------------------------------------------------------------------- <<<<<



}




void tests::MatrixXd_to_file_test(){

    unsigned int row_cnt = 10;
    unsigned int col_cnt = 10;

    Eigen::MatrixXd testMat = Eigen::MatrixXd::Random( row_cnt, col_cnt );

    string targetDir = SRC_PATH_XYQ_str + "/data_output";
    string targetStemName = "MatrixXd_to_file_res";

    utils::MatrixXd_to_file( targetDir, targetStemName, testMat, 10 );

}



void tests::file_to_MatrixXd_test(){

// ---------------------------------------------------------------------- >>>>>
//      Initialization
// ---------------------------------------------------------------------- >>>>>

    string targetDir = RES_PATH_XYQ_str + "/test_res_dir";
    string targetStemName = "test_MatrixXd_data";
    string targetExt = ".txt";
    
    Eigen::MatrixXd parsedMatAns( 10,10 );
    parsedMatAns << -9.9749748222E-01,-6.5178380688E-01,9.7705008087E-01,-6.6753135777E-01,7.5194555498E-01,5.5931272317E-01,-2.4826807459E-01,6.7522202216E-01,1.9809564501E-01,3.4211249123E-02,
    1.2717062899E-01,7.1788689840E-01,-1.0861537523E-01,3.2609027375E-01,4.5335245827E-01,6.8730735191E-01,-8.1475264748E-01,4.5298623615E-01,-2.2952970977E-01,9.7997985778E-01,
    -6.1339152196E-01,4.2100283822E-01,-7.6183355205E-01,-9.8422193060E-02,9.1180150761E-01,9.9359111301E-01,3.5441145054E-01,-3.0121768853E-02,4.7001556444E-01,5.0309762871E-01,
    6.1748100223E-01,2.7069917905E-02,-9.9066133610E-01,-2.9575487533E-01,8.5143589587E-01,9.9938962981E-01,-8.8756981109E-01,-5.8928189947E-01,2.1793267617E-01,-3.0887783441E-01,
    1.7001861629E-01,-3.9201025422E-01,-9.8217719047E-01,-8.8592181158E-01,7.8707235939E-02,2.2299874874E-01,-9.8242133854E-01,4.8747215186E-01,1.4481032746E-01,-6.6203802606E-01,
    -4.0253913999E-02,-9.7003082369E-01,-2.4423963134E-01,2.1536912137E-01,-7.1532334361E-01,-2.1512497330E-01,8.3758049257E-01,-6.3081759087E-02,-2.7732169561E-01,3.1461531419E-01,
    -2.9941709647E-01,-8.1719412824E-01,6.3325907163E-02,5.6663716544E-01,-7.5838496048E-02,-4.6757408368E-01,-4.4822534867E-01,-8.4078493606E-02,-6.9689016388E-01,-1.6205328532E-02,
    7.9192480239E-01,-2.7109591968E-01,1.4236884671E-01,6.0521256142E-01,-5.2934354686E-01,-4.0543839839E-01,-4.5420697653E-01,8.9831232643E-01,-5.4979094821E-01,-8.7292092654E-01,
    6.4568010498E-01,-7.0537430952E-01,2.0352793970E-01,3.9765617847E-02,7.2447889645E-01,6.8028809473E-01,1.7581713309E-01,4.8887600330E-01,-1.4969328898E-01,3.9951780755E-01,
    4.9320963164E-01,-6.6820276498E-01,2.1433149205E-01,-3.9609973449E-01,-5.8079775384E-01,-9.5251319926E-01,3.8236640522E-01,-7.8344065676E-01,6.0576189459E-01,9.6133304849E-03;

// ---------------------------------------------------------------------- <<<<<

    string fullFileName = targetDir + "/" + targetStemName + targetExt;
    Eigen::MatrixXd parsedMat = utils::file_to_MatrixXd( fullFileName );

    targetExt = ".csv";
    fullFileName = targetDir + "/" + targetStemName + targetExt;
    Eigen::MatrixXd parsedMat2 = utils::file_to_MatrixXd( fullFileName );

    bool test_bool = true;
    test_bool = test_bool && ( parsedMat == parsedMatAns );
    test_bool = test_bool && ( parsedMat2 == parsedMatAns );
    if( test_bool ){
        cout << "file_to_MatrixXd test: passed!" << endl;
    }else{
        cout << "file_to_MatrixXd test: failed!" << endl;
    }

}



void tests::MatrixXcd_to_file_test(){

    unsigned int row_cnt = 10;
    unsigned int col_cnt = 10;

    Eigen::MatrixXd realPart = utils::gen_rand_MatrixXd( row_cnt, col_cnt );
    Eigen::MatrixXd imagPart = utils::gen_rand_MatrixXd( row_cnt, col_cnt );

    Eigen::MatrixXcd testMat( row_cnt, col_cnt );
    testMat.real() = realPart;
    testMat.imag() = imagPart;

    string targetDir = SRC_PATH_XYQ_str + "/data_output";
    string targetStemName = "MatrixXcd_to_file_res";


    utils::MatrixXcd_to_file( targetDir, targetStemName, testMat, 10 );

}
