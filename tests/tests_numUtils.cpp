#include "tests_numUtils.h"


using namespace std;


void tests::numUtils_test_1( unsigned int case_idx ){

    // Initialize test case index.
    int case_cnt = 0;

    // 0- Test gen_lin_idx_arr.
    if( case_cnt == case_idx ){

        bool match_flag = true;

        vector< unsigned int > testVec = utils::gen_lin_idx_arr( 0, 10, 11 );
        vector< unsigned int > ansVec = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        if( testVec.size() == ansVec.size() ){
            for( unsigned int z = 0; z < ansVec.size(); z++ ){
                match_flag = match_flag && ( testVec.at(z) == ansVec.at(z) );
            }
        }else{
            match_flag = false;
        }

        cout << "Check 1 (Full index case): " << match_flag << endl;

    }


    case_cnt++;
    // 1- Test gen_rem_idx_arr.
    if( case_cnt == case_idx ){
        
        bool match_flag = true;

        vector< unsigned int > p1vec = utils::gen_lin_idx_arr( 1, 17, 3 );
        vector< unsigned int > p2vec = utils::gen_rem_idx_arr( 0, 20, p1vec );
        vector< unsigned int > p2vec_ans = { 0, 2, 3, 4, 5, 6, 7, 8, 10, 11, 
            12, 13, 14, 15, 16, 18, 19, 20 };
        
        if( p2vec.size() == p2vec_ans.size() ){
            for( unsigned int z = 0; z < p2vec_ans.size(); z++ ){
                match_flag = match_flag && ( p2vec.at(z) == p2vec_ans.at(z) );
            }
        }else{
            match_flag = false;
        }
        cout << "Check 1: " << match_flag << endl;

    }

    case_cnt++;
    // 2- Test gen_even_idx_arr and gen_odd_idx_arr.
    if( case_cnt == case_idx ){
        
        bool match_flag = true;

        vector< unsigned int > p1vec = utils::gen_even_idx_arr( 1, 17 );
        vector< unsigned int > p2vec = utils::gen_odd_idx_arr( 1, 17 );
        vector< unsigned int > p1vec_ans = { 2, 4, 6, 8, 10, 12, 14, 16 };
        vector< unsigned int > p2vec_ans = { 1, 3, 5, 7, 9, 11, 13, 15, 17};

        if( p1vec.size() == p1vec_ans.size() ){
            for( unsigned int z = 0; z < p1vec_ans.size(); z++ ){
                match_flag = match_flag && ( p1vec.at(z) == p1vec_ans.at(z) );
            }
        }else{
            match_flag = false;
        }
        cout << "Check even array: " << match_flag << endl;

        match_flag = true;
        if( p2vec.size() == p2vec_ans.size() ){
            for( unsigned int z = 0; z < p2vec_ans.size(); z++ ){
                match_flag = match_flag && ( p2vec.at(z) == p2vec_ans.at(z) );
            }
        }else{
            match_flag = false;
        }
        cout << "Check odd array: " << match_flag << endl;

    }


    case_cnt++;
    // 3- Test rIntGen()
    if( case_cnt == case_idx ){

        bool match_flag = true;
        int L_bnd = 3;
        int U_bnd = 26;
        unsigned int cnt = 10;

        // Generate a random integer array.
        shared_ptr<vector<int>> myRandIntVec = utils::rIntGen( L_bnd, U_bnd, cnt );

        // Find the maximum element using std::max_element
        auto maxElemIter = std::max_element(myRandIntVec->begin(), myRandIntVec->end());
        auto minElemIter = std::min_element(myRandIntVec->begin(), myRandIntVec->end());

        match_flag = match_flag && ( myRandIntVec->size() == cnt );
        for( unsigned int z = 0; z < cnt; z++ ){
            match_flag = match_flag && ( myRandIntVec->at(z) <= *maxElemIter );
            match_flag = match_flag && ( myRandIntVec->at(z) >= *minElemIter );
            if( !match_flag ){
                break;
            }
        }

        cout << "Random integer generator test match: " << match_flag << endl;

    }


    case_cnt++;
    // 4- Test rDoubleGen()
    if( case_cnt == case_idx ){

        bool match_flag = true;
        double L_bnd = 1.245;
        double U_bnd = 67.293;
        unsigned int cnt = 20;

        // Generate a random integer array.
        shared_ptr<vector<double>> myRandIntVec = utils::rDoubleGen( L_bnd, U_bnd, cnt );

        // Find the maximum element using std::max_element
        auto maxElemIter = std::max_element(myRandIntVec->begin(), myRandIntVec->end());
        auto minElemIter = std::min_element(myRandIntVec->begin(), myRandIntVec->end());

        match_flag = match_flag && ( myRandIntVec->size() == cnt );
        for( unsigned int z = 0; z < cnt; z++ ){
            match_flag = match_flag && ( myRandIntVec->at(z) <= *maxElemIter );
            match_flag = match_flag && ( myRandIntVec->at(z) >= *minElemIter );
            if( !match_flag ){
                break;
            }
        }

        cout << "Random double generator test match: " << match_flag << endl;

    }


}


void tests::gen_match_vector_test(){


    vector<unsigned int> vec_1A = { 0, 1, 2, 6, 7, 7, 3 };
    vector<unsigned int> vec_1B = { 1, 2, 7, 11, 5, 7 };
    vector<unsigned int> vec_1C = utils::gen_match_vector( vec_1A, vec_1B );

    vector<unsigned int> vec_1C_ans = { 1, 2, 7, 7 };
    bool test_bool = true;
    test_bool = test_bool && ( vec_1C_ans.size() == vec_1C.size() );
    if( test_bool ){
        for( unsigned int z = 0; z < vec_1C_ans.size(); z++ ){
            test_bool = test_bool && ( vec_1C_ans[z] == vec_1C[z] );
        }
    }
    if( test_bool ){
        cout << "gen_match_vector_test test 1: passed!" << endl;
    }else{
        cout << "gen_match_vector_test test 1: failed!" << endl;
    }


    vector<unsigned int> vec_2A = { 0, 1, 2, 4 };
    vector<unsigned int> vec_2B = { 6, 7, 9 };
    vector<unsigned int> vec_2C = utils::gen_match_vector( vec_2A, vec_2B );

    test_bool = true;
    test_bool = test_bool && ( vec_2C.size() == 0 );
    if( test_bool ){
        cout << "gen_match_vector_test test 2: passed!" << endl;
    }else{
        cout << "gen_match_vector_test test 2: failed!" << endl;
    }



}


void tests::sort_num_vec_inplace_test(){

    vector<unsigned int> vec_1A_ans = { 0, 1, 2, 3, 5, 8, 9 };
    vector<unsigned int> vec_1A = { 0, 1, 8, 3, 2, 9, 5 };
    utils::sort_num_vec_inplace( vec_1A, true );
    bool test_bool = true;
    if( test_bool ){
        for( unsigned int z = 0; z < vec_1A.size(); z++ ){
            test_bool = test_bool && ( vec_1A[z] == vec_1A_ans[z] );
        }
    }
    if( test_bool ){
        cout << "gen_match_vector_test test 1: passed!" << endl;
    }else{
        cout << "gen_match_vector_test test 1: failed!" << endl;
    }

    vector<double> vec_2A_ans = { 0.982, 0.98, 0.75, 0.71, 0.566, 0.43, 0.2 };
    vector<double> vec_2A =     { 0.43, 0.2, 0.75, 0.982, 0.71, 0.566, 0.98 };
    utils::sort_num_vec_inplace( vec_2A, false );
    test_bool = true;
    if( test_bool ){
        for( unsigned int z = 0; z < vec_2A.size(); z++ ){
            test_bool = test_bool && ( abs( vec_2A[z] - vec_2A_ans[z] ) < 1e-12 );
        }
    }
    if( test_bool ){
        cout << "gen_match_vector_test test 2: passed!" << endl;
    }else{
        cout << "gen_match_vector_test test 2: failed!" << endl;
    }

}


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

