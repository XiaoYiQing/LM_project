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

}

