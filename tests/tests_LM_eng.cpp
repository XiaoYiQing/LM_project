#include "tests_LM_eng.h"



using namespace std;

void tests::LM_eng_test_1( unsigned int test_idx ){


    // Define our frequency data object.
    fData myF;

    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
    fData::read_sXp_file( myF, fullFileName );

    vector< unsigned int > fr_idx_arr = utils::gen_lin_idx_arr( 0, myF.get_f_cnt() - 1, 10 );
    shared_ptr<fData> myFr = myF.red_partit( fr_idx_arr );

    // Generate two partitions based on linear distribution.
    vector< shared_ptr<fData> > myPartits = myFr->gen_2_partit();

    shared_ptr<fData> partit1 = myPartits.at(0);
    shared_ptr<fData> partit2 = myPartits.at(1);

    Eigen::MatrixXcd myLM = LM_UTIL::build_LM( *partit1, *partit2 );

}
