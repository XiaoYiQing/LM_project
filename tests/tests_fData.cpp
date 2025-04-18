#include "tests_fData.h"


using namespace std;

void tests::fData_test_1( unsigned int test_idx ){


    int case_cnt = 0;


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


}
