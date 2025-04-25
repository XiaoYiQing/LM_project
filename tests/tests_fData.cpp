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

        


    }


}
