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


void tests::file_to_MatrixXcd_test(){

// ---------------------------------------------------------------------- >>>>>
//      Initialization
// ---------------------------------------------------------------------- >>>>>

    string targetDir = RES_PATH_XYQ_str + "/test_res_dir";
    string targetStemName = "test_MatrixXcd_data";
    string targetExt = ".txt";

    Eigen::MatrixXd matAns_real( 10,10 );
    matAns_real << 9.2228925733E-01,	6.9808128987E-01,	4.6199454121E-01,	3.0746261246E-01,	-5.4423972415E-01,	5.6793968396E-01,	-8.9463791910E-01,	-3.9603774225E-01,	-1.2061888943E-01,	6.3805883082E-01,
    6.1064886408E-01,	-1.5556741450E-01,	5.9629913990E-01,	8.8437924663E-01,	-1.5152034074E-01,	1.9577702023E-03,	-2.2213529529E-02,	-7.3448016562E-01,	2.2647690107E-01,	-6.2574193966E-01,
    -4.5554897561E-01,	-8.5396212442E-01,	-9.0578575166E-01,	-8.1916771358E-02,	7.7633702912E-01,	4.0945773338E-01,	-7.8183709294E-01,	-9.3397789983E-01,	-2.9750373272E-01,	7.1874041857E-01,
    -6.8728193350E-01,	3.1544475112E-01,	4.4298740785E-01,	-8.5588957444E-01,	6.9048850831E-01,	-7.8469591364E-02,	1.4397391402E-01,	8.6011714930E-01,	-2.3563231825E-02,	-5.2550682399E-01,
    -3.7011660266E-01,	3.7779981860E-01,	-4.6491719856E-01,	-2.0022581508E-01,	5.0727558592E-02,	6.4926756551E-01,	8.0030173968E-01,	-5.8640725270E-01,	-1.6012881877E-01,	-6.2260601392E-01,
    -4.6497094039E-01,	-8.4671119479E-01,	-7.1347666715E-01,	-9.2791787700E-01,	9.7983504166E-01,	9.1811035230E-01,	3.1417221561E-01,	5.3749247819E-01,	3.0494587086E-01,	-8.3840932435E-01,
    -7.3291448474E-01,	7.8706672449E-01,	-3.3475658910E-01,	8.8310469967E-01,	-2.5964243377E-01,	-8.8136728178E-01,	-6.2724950058E-01,	-4.8147101223E-01,	3.0778289859E-01,	-4.5144475226E-01,
    9.3866751128E-01,	1.5126364311E-01,	5.2973109086E-01,	-2.1590733506E-01,	-9.0155733852E-01,	-6.0621657247E-02,	3.4989936923E-01,	2.7165527026E-02,	-4.4787178389E-01,	6.0525769256E-01,
    9.1408008960E-01,	-5.1883734697E-02,	-6.3682573697E-01,	-5.2799871442E-01,	4.7255016763E-01,	8.4555062550E-01,	-2.4716006873E-02,	4.9787759898E-01,	4.1364957967E-01,	-2.4278037877E-01,
    -5.1964916733E-01,	-3.5776164397E-01,	9.8966002943E-01,	-6.5532639655E-01,	8.5865151784E-01,	3.6623312942E-01,	-6.7728167953E-01,	-1.4158185567E-01,	-4.8520294953E-01,	-1.5161479035E-01;

    Eigen::MatrixXd matAns_imag( 10,10 );
    matAns_imag << -6.5217710889E-01,	1.9099359419E-01,	-2.6110970487E-01,	-7.2128480502E-03,	-4.9932012491E-01,	4.9973369641E-01,	-2.1654363177E-01,	-8.5831323149E-01,	-6.3045149432E-01,	-6.9783500334E-01,
    7.3467543381E-01,	-3.2250757228E-01,	-2.8309708043E-01,	7.5989011207E-01,	-2.2033649116E-01,	2.9852523951E-01,	1.2065205755E-01,	2.1276250242E-01,	-1.5473291762E-01,	4.8158494047E-01,
    -6.2646741610E-01,	1.8735730769E-01,	8.3445211482E-01,	-2.7100106170E-02,	8.0023287448E-01,	9.4871991896E-01,	7.7823728150E-01,	-4.9159842487E-01,	-8.1799584570E-01,	-6.6729961318E-03,
    9.7637580601E-01,	-2.6053371621E-01,	7.5521221862E-01,	-2.5685755974E-02,	4.5063704760E-01,	7.2157685568E-01,	9.6932130069E-01,	-6.6989637337E-01,	2.7573603519E-01,	6.3267892325E-02,
    -6.5709341607E-01,	-4.2978731439E-01,	-5.1375776357E-02,	5.2119518188E-01,	-3.7126127723E-01,	8.6856225792E-01,	-8.7600489022E-01,	-6.4642617287E-02,	7.0521190991E-01,	-1.2780870586E-01,
    8.8697963370E-01,	1.7499686591E-01,	-8.5732291534E-02,	-4.1835268691E-01,	2.2655135764E-01,	-4.1392211165E-01,	6.1181440042E-02,	-1.4526024445E-01,	6.7616983583E-01,	9.7759889816E-01,
    7.1097965853E-01,	-9.8018082471E-01,	-3.5418613584E-01,	8.5592857849E-01,	2.7034019745E-01,	8.3913878275E-01,	-5.4178861788E-01,	9.6922168699E-01,	-6.1669626846E-01,	1.4174824366E-01,
    8.2069491247E-01,	-3.4658502792E-01,	-5.9205164584E-01,	1.0139971217E-01,	5.9611661145E-01,	-4.6206396936E-01,	6.1281026577E-01,	-7.4194717361E-01,	7.3979076401E-01,	-6.2086008003E-01,
    -9.4216805162E-01,	9.3393994059E-01,	-8.4088653071E-01,	5.5268130012E-01,	-7.2234635596E-01,	4.9778063212E-02,	-5.6022231954E-01,	-5.2183197470E-01,	1.5137262866E-01,	1.2482308439E-01,
    -2.8723276180E-01,	4.6923930302E-01,	8.9785385712E-01,	-1.1931896938E-01,	4.4952273371E-01,	9.0182976852E-01,	6.5884896247E-01,	-1.8948503885E-01,	3.6278199911E-02,	4.2264630488E-01;

    Eigen::MatrixXcd matAns( 10, 10 );
    matAns.real() = matAns_real;
    matAns.imag() = matAns_imag;

// ---------------------------------------------------------------------- <<<<<

    string fullFileName = targetDir + "/" + targetStemName + targetExt;
    Eigen::MatrixXcd parsedMat = utils::file_to_MatrixXcd( fullFileName );

    targetExt = ".csv";
    fullFileName = targetDir + "/" + targetStemName + targetExt;
    Eigen::MatrixXcd parsedMat2 = utils::file_to_MatrixXcd( fullFileName );

    bool test_bool = true;
    test_bool = test_bool && ( parsedMat == matAns );
    test_bool = test_bool && ( parsedMat2 == matAns );
    if( test_bool ){
        cout << "file_to_MatrixXcd test: passed!" << endl;
    }else{
        cout << "file_to_MatrixXcd test: failed!" << endl;
    }

}



void tests::gen_rand_MatrixXd_test(){

    Eigen::MatrixXd tmp = utils::gen_rand_MatrixXd( 4, 4 );
    cout << tmp << endl;

}


void tests::gen_randn_MatrixXd_test(){

    bool test_bool = true;

    unsigned int row_cnt = 1000;
    unsigned int col_cnt = 100;
    Eigen::MatrixXd tmp = utils::gen_randn_MatrixXd( row_cnt, col_cnt );

    // Transfer all values of the matrix into a vector.
    Eigen::VectorXd entriesVec( row_cnt * col_cnt );
    for( unsigned int z = 0; z < row_cnt*col_cnt; z++ ){
        entriesVec(z) = tmp.data()[z];
    }

    // Compute mean.
    double mean = entriesVec.mean();
    // Calculate variance.
    double variance = (entriesVec.array() - mean).square().mean();
    // Compute standard deviation.
    double std_dev = std::sqrt(variance);

    cout << "Mean (Expect 0): " << mean << endl;
    cout << "Standard deviation (Expect 1): " << std_dev << endl;

}


void tests::gen_orth_basis_test(){

    // Define the numerical threshold.
    double num_thresh = 1e-12;

    // Define the dimnesions of the test matrix.
    unsigned int row_cnt = 10;
    unsigned int col_cnt = 10;

    // Define a random test matrix.
    Eigen::MatrixXd testMat = utils::gen_rand_MatrixXd( row_cnt, col_cnt );
    // Generate an orthonormal basis.
    Eigen::MatrixXd Q_mat = utils::gen_orth_basis( testMat );

    // Initialize test boolean.
    bool test_bool = true;
    // Initialize two vectors.
    Eigen::VectorXd v1;
    Eigen::VectorXd v2;

    for( unsigned int i = 0; i < row_cnt; i++ ){
        v1 = Q_mat.col(i);
        for( unsigned int j = 0; j < col_cnt; j++ ){
            v2 = Q_mat.col(j);

            // Dot product verification for orthonormality.
            if( i == j ){
                test_bool = test_bool && ( abs( v1.dot(v2) - 1 ) < num_thresh );
            }else{
                test_bool = test_bool && ( abs( v1.dot(v2) ) < num_thresh );
            }

        }
    }

    if( test_bool ){
        cout << "gen_orth_basis test: passed!" << endl;
    }else{
        cout << "gen_orth_basis test: failed!" << endl;
    }

}


void tests::SVD_econ_test(){

    // Define the numerical threshold.
    double num_thresh = 1e-12;

    // Define the dimensions of the test matrix.
    unsigned int row_cnt = 300;
    unsigned int col_cnt = 30;


    // Define a random test matrix.
    auto start = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXd testMat = utils::gen_rand_MatrixXd( row_cnt, col_cnt );
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Econ SVD runtime: " << duration.count() << " ms" << std::endl;
    

    // Initialize test boolean.
    bool test_bool = true;

    // Perform the 'econ' SVD.
    utils::SVD_econ SVD_res( testMat );
    // Obtain the results of the SVD.
    Eigen::MatrixXd U = SVD_res.U;
    Eigen::MatrixXd V = SVD_res.V;
    Eigen::VectorXd S = SVD_res.S;
    // Check key property of the singular vectors.
    Eigen::MatrixXd Ut_x_U_res = U.transpose()*U;
    Eigen::MatrixXd Vt_x_V_res = V.transpose()*V;

    // Singular vectors multiplication results check.
    test_bool = test_bool && ( U.rows() == row_cnt && U.cols() == col_cnt );
    test_bool = test_bool && ( V.rows() == col_cnt && V.cols() == col_cnt );
    test_bool = test_bool && ( S.size() == col_cnt );
    for( unsigned int i = 0; i < col_cnt; i++ ){
    for( unsigned int j = 0; j < col_cnt; j++ ){
        if( i == j ){
            test_bool = test_bool && ( abs( Ut_x_U_res(i,j) - 1 ) < num_thresh );
            test_bool = test_bool && ( abs( Vt_x_V_res(i,j) - 1 ) < num_thresh );
        }else{
            test_bool = test_bool && ( abs( Ut_x_U_res(i,j) ) < num_thresh );
            test_bool = test_bool && ( abs( Vt_x_V_res(i,j) ) < num_thresh );
        }
    }
    }

    if( test_bool ){
        cout << "SVD_econ test: passed!" << endl;
    }else{
        cout << "SVD_econ test: failed!" << endl;
    }

    auto start2 = std::chrono::high_resolution_clock::now();
    Eigen::JacobiSVD<Eigen::MatrixXd> full_SVD( testMat, Eigen::ComputeFullU | Eigen::ComputeFullV );
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration2 = end2 - start2;
    std::cout << "Full SVD runtime: " << duration2.count() << " ms" << std::endl;

}
