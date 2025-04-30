#include "tests_LM_eng.h"



using namespace std;

void tests::LM_eng_test_1( unsigned int test_idx ){


    // Define our frequency data object.
    fData myF;

    // Define the full file name.
    string fullFileName = RES_PATH_XYQ_str + "/Slink_a=100um_b=400um.s2p";
    fData::read_sXp_file( myF, fullFileName );

    // Switch the data format into real + imaginary format.
    myF.data_format_Switch( fData::FDATA_FORMAT::RI );

    // Create a subset linear index array.
    vector< unsigned int > fr_idx_arr = utils::gen_lin_idx_arr( 0, myF.get_f_cnt() - 1, 100 );
    // Create a fData subset.
    shared_ptr<fData> myFr = myF.red_partit( fr_idx_arr );

    // Generate two partitions from this data subset.
    vector< shared_ptr<fData> > myPartits = myFr->gen_2_partit();


    shared_ptr<fData> partit1 = myPartits.at(0);
    shared_ptr<fData> partit2 = myPartits.at(1);
    // Constrcut the Loewner Matrix using the two partitions.
    Eigen::MatrixXcd myLM = LM_UTIL::build_LM( *partit1, *partit2 );


    // Compute the SVD
    Eigen::JacobiSVD<Eigen::MatrixXcd> mySVD(myLM, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Get the singular values
    Eigen::VectorXd singularValues = mySVD.singularValues();

    cout << singularValues << endl;

    int lol = 0;

}
