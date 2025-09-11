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

        string targetDir = "C:/Users/Yi Qing Xiao/Documents/Cpp_projects/LM_project/data_output";
        string targetStemName = "tmp_data_file";

        utils::vec_to_file( targetDir, targetStemName, tarVec, 0 );

    }

}


void tests::file_to_vec_test(){

    string targetDir = RES_PATH_XYQ_str + "/test_res_dir";
    string targetStemName = "test_vec_data";
    string targetExt = ".txt";

    string fullFileName = targetDir + "/" + targetStemName + targetExt;

    vector<double> my_vec = utils::file_to_vec( fullFileName );
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
        cout << "file_to_vec correct read test: passed!" << endl;
    }else{
        cout << "file_to_vec correct read test: failed!" << endl;
    }
    cout << endl;

    targetExt = ".csv";
    fullFileName = targetDir + "/" + targetStemName + targetExt;

    vector<double> my_vec2 = utils::file_to_vec( fullFileName );

    test_vect = my_vec2.size() == my_vec_ans.size();
    for( unsigned int z = 0; z < my_vec2.size(); z++ ){
        test_vect = test_vect && ( abs( my_vec2[z] - my_vec_ans[z] ) < 1e-12 );
    }
    if( test_vect ){
        cout << "file_to_vec correct csv file read test: passed!" << endl;
    }else{
        cout << "file_to_vec correct csv file read test: failed!" << endl;
    }
    cout << endl;




    targetDir = RES_PATH_XYQ_str + "/test_res_dir";
    targetStemName = "test_bad_vec_data";
    targetExt = ".txt";

    fullFileName = targetDir + "/" + targetStemName + targetExt;
    try{
        vector<double> my_vec3 = utils::file_to_vec( fullFileName );
        test_vect = false;
    }catch(...){
        test_vect = true;
    }
    if( test_vect ){
        cout << "file_to_vec bad file read test: passed!" << endl;
    }else{
        cout << "file_to_vec bad file read test: failed!" << endl;
    }
    cout << endl;

}
