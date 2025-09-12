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

        utils::VectorXd_to_file( targetDir, targetStemName, tarVec, 0 );

    }

}


void tests::Vd_to_file_test(){

    unsigned int num_cnt = 10;

    vector<double> tmp_vec = *utils::rDoubleGen( -10.0, 10.0, num_cnt );

    string targetDir = SRC_PATH_XYQ_str + "/data_output";
    string targetStemName = "Vd_to_file_res";

    utils::Vd_to_file( targetDir, targetStemName, tmp_vec, 0 );
    
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

    utils::MatrixXd_to_file( targetDir, targetStemName, testMat, 0 );

}



void tests::file_to_MatrixXd(){

// ---------------------------------------------------------------------- >>>>>
//      Initialization
// ---------------------------------------------------------------------- >>>>>

    string targetDir = RES_PATH_XYQ_str + "/test_res_dir";
    string targetStemName = "test_MatrixXd_data";
    string targetExt = ".txt";
    string fullFileName = targetDir + "/" + targetStemName + targetExt;

// ---------------------------------------------------------------------- <<<<<

    utils::file_to_MatrixXd( fullFileName );

}
