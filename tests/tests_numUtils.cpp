#include "tests_numUtils.h"


using namespace std;


void tests::numUtils_test( unsigned int case_idx ){

    // Initialize test case index.
    int case_cnt = 0;

    // 0- Empty vector initialization exception case.
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

}

