#include "tests_Matrix3DXd.h"

using namespace std;


void tests::Matrix3DXd_test_1( int case_idx ){

    // Initialize test case index.
    int case_cnt = 0;


    // 0- Empty vector initialization exception case.
    if( case_cnt == case_idx ){

        vector< Eigen::MatrixXd > tmp_mat_vec;

        try{
            Matrix3DXd my_3D_mat = Matrix3DXd( tmp_mat_vec );
            cout << "Empty vector initialization test: Failed!" << endl;
        }catch ( const std::out_of_range& e ){
            cerr << e.what() << endl;
            cout << "Empty vector initialization test: Passed!" << endl;
            // return;
        }

        Matrix3DXd my_3D_mat2 = Matrix3DXd( 0, 0, 0 );
        bool test_1_check = my_3D_mat2.isEmpty() == true;
        test_1_check = test_1_check && my_3D_mat2.rows() == 0;
        test_1_check = test_1_check && my_3D_mat2.cols() == 0;
        test_1_check = test_1_check && my_3D_mat2.levels() == 0;
        cout << "Emtpy initialization test: ";
        if( test_1_check ){
            cout << "Passed!" << endl;
        }else{
            cout << "Failed!" << endl;
        }

    }

    case_cnt++;
    // 1- Multiplication operator test.
    if( case_cnt == case_idx ){

        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }

        unsigned int test_idx = 2;
        double multiplier = 3.0;
        Matrix3DXd my3DMat = Matrix3DXd( tmp_mat_vec );
        Matrix3DXd my3DMat2 = my3DMat*multiplier;

        cout << "Matrix3DXd scalar multiplcation operator test: ";
        if(  multiplier*tmp_mat_vec.at( test_idx ) == my3DMat2.at( test_idx ) ){
            cout << "Passed!" << endl;
        }else{
            cout << "Failed!" << endl;
        }

                
    }
    
    case_cnt++;
    // 2- Inplace multiplication operator test.
    if( case_cnt == case_idx ){

        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z+0.1,z+0.2,z+0.3,z+0.4;
            tmp_mat_vec.push_back( tmpMat );
        }

        Matrix3DXd my3DMat = Matrix3DXd( tmp_mat_vec );

        bool test_bool = true;
        double factor = 3;
        my3DMat *= factor;
        for( unsigned int z = 0; z < tmp_mat_vec.size(); z++ ){

            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z+0.1,z+0.2,z+0.3,z+0.4;
            tmpMat *= factor;

            test_bool = test_bool && ( my3DMat.at(z) == tmpMat );

        }

        if( test_bool ){
            cout << "Matrix3DXd inplace scalar multiplcation operator test: passed!" << endl;
        }else{
            cout << "Matrix3DXd inplace scalar multiplcation operator test: failed!" << endl;
        }
        
                
    }


    case_cnt++;
    // 3- Invalid input vector initialization check.
    if( case_cnt == case_idx ){

        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }
        // Add another matrix to the vector, but having different dimensions 
        // than previous matrices.
        Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 3 );
        tmpMat << 9,9,9,9,9,9;
        tmp_mat_vec.push_back( tmpMat );

        bool test_bool = true;
        try{
            Matrix3DXd my_3D_mat = Matrix3DXd( tmp_mat_vec );
            test_bool = false;
        }catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
        }
        if( test_bool ){
            cout << "Matrix3DXd push_back wrong sized entry test: passed!" << endl;
        }else{
            cout << "Matrix3DXd push_back wrong sized entry test: failed!" << endl;
        }

    }


    case_cnt++;
    // 4- push_back function check.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }
        
        // Use the vector of 2D matrices to initialize a 3D matrix.
        Matrix3DXd my_3D_mat;
        try{
            my_3D_mat = Matrix3DXd( tmp_mat_vec );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        }
        

        bool test_bool = true;

        // Add another matrix to the vector at its end.
        Eigen::MatrixXd addMat = Eigen::MatrixXd( 2, 2 );
        addMat << 91, 92, 93, 94;
        try{
            my_3D_mat.push_back( addMat );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            test_bool = false;
        }
        test_bool = test_bool && ( addMat == my_3D_mat.at( my_3D_mat.levels() - 1 ) );

        if( test_bool ){
            cout << "Matrix3DXd push_back correct sized entry test: passed!" << endl;
        }else{
            cout << "Matrix3DXd push_back correct sized entry test: failed!" << endl;
        }

    }



    case_cnt++;
    // 5- insert function check.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }

        // Use the vector of 2D matrices to initialize a 3D matrix.
        Matrix3DXd my_3D_mat;
        try{
            my_3D_mat = Matrix3DXd( tmp_mat_vec );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        }
        
        // Insert another matrix to the vector at the target valid index.
        Eigen::MatrixXd addMat = Eigen::MatrixXd( 2, 2 );
        addMat << 91, 92, 93, 94;
        unsigned int tarIdx = 2;
        try{
            my_3D_mat.insert( tarIdx, addMat );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        } catch( const std::out_of_range& e ){
            cerr << e.what() << endl;
            return;
        }
        bool test_bool = true;
        test_bool = test_bool && ( my_3D_mat.at( tarIdx ) == addMat );
        if( test_bool ){
            cout << "Matrix3DXd insert correct sized entry test: passed!" << endl;
        }else{
            cout << "Matrix3DXd insert correct sized entry test: failed!" << endl;
        }


        // Insert another matrix to the vector at the target invalid index.
        Eigen::MatrixXd addMat2 = Eigen::MatrixXd( 2, 2 );
        addMat2 << 51, 52, 53, 54;
        tarIdx = 100;
        try{
            my_3D_mat.insert( tarIdx, addMat2 );
            test_bool = false;
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            test_bool = true;
        } catch( const std::out_of_range& e ){
            cerr << e.what() << endl;
            test_bool = true;
        } catch( ... ){
            test_bool = false;
        }

        if( test_bool ){
            cout << "Matrix3DXd insert incorrect index test: passed!" << endl;
        }else{
            cout << "Matrix3DXd insert incorrect index test: failed!" << endl;
        }

    }


    case_cnt++;
    // 6- Batch push_back function check.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }
        // Create the 3D matrix object using the above matrix vector.
        Matrix3DXd my_3D_mat = Matrix3DXd( tmp_mat_vec );

        vector< Eigen::MatrixXd > addVec;
        unsigned int add_cnt = 3;
        for( unsigned int z = 0; z < add_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z*10+1,z*10+2,z*10+3,z*10+4;
            addVec.push_back( tmpMat );
        }

        // Put the additional matrix vector at the end of the 3D matrix.
        my_3D_mat.push_back( addVec );
        cout << my_3D_mat.at( my_3D_mat.levels() - 1 ) << endl;

    }


    case_cnt++;
    // 7- Reserve function check.
    if( case_cnt == case_idx ){

        // Create the 3D matrix object using the above matrix vector.
        Matrix3DXd my_3D_mat = Matrix3DXd();


        Eigen::MatrixXd init_mat = Eigen::MatrixXd( 2, 2 );
        init_mat << 1, 2, 3, 4;
        my_3D_mat.push_back( init_mat );

        cout << "3D Matrix levels before reserve: " << my_3D_mat.levels() << endl;
        my_3D_mat.reserve( 200 );
        cout << "3D Matrix levels after reserve: " << my_3D_mat.levels() << endl;


    }

    case_cnt++;
    // 8- set() and at() functions check.
    if( case_cnt == case_idx ){
        // Create the 3D matrix object using the above matrix vector.
        Matrix3DXd my_3D_mat = Matrix3DXd();

        Eigen::MatrixXd init_mat = Eigen::MatrixXd( 2, 2 );
        init_mat << 1, 2, 3, 4;
        my_3D_mat.push_back( init_mat );

        Eigen::MatrixXd next_mat = Eigen::MatrixXd( 2, 2 );
        next_mat << 11, 12, 13, 14;
        my_3D_mat.at(0) = next_mat;
        cout << my_3D_mat.at(0) << endl;
        my_3D_mat.set( 0, next_mat );
        cout << my_3D_mat.at(0) << endl;
    }


    case_cnt++;
    // 9- at( vector< unsigned int > idxVec ).
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 9;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 10*z+1,10*z+2,10*z+3,10*z+4;
            tmp_mat_vec.push_back( tmpMat );
        }
        // Initialize a 3D matrix using the vector.
        Matrix3DXd my3DMat = Matrix3DXd( tmp_mat_vec );


        vector<unsigned int> sub_vec_idx = { 2, 3, 4 }; 
        Matrix3DXd my3DMatSub = my3DMat.at( sub_vec_idx );
        bool match = true;
        for( unsigned int z = 0; z < sub_vec_idx.size(); z++ ){

            unsigned int curr_idx = sub_vec_idx[z];

            match = match && ( my3DMatSub.at(z) == my3DMat.at( curr_idx ) );

        }

        cout << "Matching data? -> " << match << endl;

    }

    return;


}




void tests::Matrix3DXd_test_2_ops( int case_idx ){

    // Initialize test case index.
    int case_cnt = 0;

    // 0- elem_pow test.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 10*z+1,10*z+2,10*z+3,10*z+4;
            tmp_mat_vec.push_back( tmpMat );
        }

        // Initialize a 3D matrix using the vector.
        Matrix3DXd my3DMat = Matrix3DXd( tmp_mat_vec );

        cout << my3DMat.at( 3 ) << endl;
        my3DMat.elem_pow( 3 );
        cout << my3DMat.at( 3 ) << endl;
        my3DMat.elem_pow( 1.0/3 );
        cout << my3DMat.at( 3 ) << endl;

    }


    case_cnt++;
    // 1- elem_raise_pow test.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 0.1+z,0.2+z,0.3+z,0.4+z;
            tmp_mat_vec.push_back( tmpMat );
        }

        // Initialize a 3D matrix using the vector.
        Matrix3DXd my3DMat = Matrix3DXd( tmp_mat_vec );

        cout << my3DMat.at( 3 ) << endl;
        my3DMat.elem_raise_pow( 2 );
        cout << my3DMat.at( 3 ) << endl;

    }

    case_cnt++;
    // 2- Same type mat element-wise mult. operator test.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 1+10*z, 2+10*z, 3+10*z, 4+10*z;
            tmp_mat_vec.push_back( tmpMat );
        }
        // Initialize a 3D matrix using the vector.
        Matrix3DXd my3DMat1 = Matrix3DXd( tmp_mat_vec );

        vector< Eigen::MatrixXd > tmp_mat_vec2;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 1.5+z, 2.5+z, 3.5+z, 4.5+z;
            tmp_mat_vec2.push_back( tmpMat );
        }
        // Initialize a 3D matrix using the vector.
        Matrix3DXd my3DMat2 = Matrix3DXd( tmp_mat_vec2 );

        // Perform matrix multiplication.
        Matrix3DXd my3DMat_prod = my3DMat1 * my3DMat2;
        cout << my3DMat1.at(1) << endl;
        cout << my3DMat2.at(1) << endl;
        cout << my3DMat_prod.at(1) << endl;

    }

    case_cnt++;
    // 3- elem_div_spec check.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 1+10*z, 2+10*z, 3+10*z, 4+10*z;
            tmp_mat_vec.push_back( tmpMat );
        }
        // Initialize a 3D matrix using the vector.
        Matrix3DXd my3DMat1 = Matrix3DXd( tmp_mat_vec );

        vector< Eigen::MatrixXd > tmp_mat_vec2;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 1.5+z, 2.5+z, 3.5+z, 4.5+z;
            tmp_mat_vec2.push_back( tmpMat );
        }
        // Initialize a 3D matrix using the vector.
        Matrix3DXd my3DMat2 = Matrix3DXd( tmp_mat_vec2 );

        // Perform element-wise matrix division.
        cout << my3DMat1.at(1) << endl;
        cout << my3DMat2.at(1) << endl;
        Matrix3DXd my3DMat_quot = my3DMat1 / my3DMat2;
        cout << my3DMat_quot.at(1) << endl;


        // Perform modification on the divisor's entry to have 0.
        my3DMat2.set( 0, 1, 1, 0 );

        Matrix3DXd resMat;

        cout << endl;
        // Perform sepcial element-wise matrix division.
        cout << my3DMat1.at(1) << endl;
        cout << my3DMat2.at(1) << endl;
        try{
            resMat = my3DMat1.elem_div_spec( my3DMat2 );
            cout << resMat.at(1) << endl;
        }catch( std::runtime_error e ){
            cerr << e.what() << endl;
        }
        

        // Perform another modification to put the same coordiante target
        // in the dividend to zero.
        my3DMat1.set( 0, 1, 1, 0 );
        cout << endl;
        // Perform sepcial element-wise matrix division.
        cout << my3DMat1.at(1) << endl;
        cout << my3DMat2.at(1) << endl;
        try{
            resMat = my3DMat1.elem_div_spec( my3DMat2 );
            cout << resMat.at(1) << endl;
        }catch( std::runtime_error e ){
            cerr << e.what() << endl;
        }
        

    }

}



void tests::Matrix3DXd_test_3_spec_ops( int case_idx ){

    // Initialize test case index.
    int case_cnt = 0;

    // 0- RMS_total_comp test.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 2;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 0+z, 1+z, 2+z, 3+z;
            tmp_mat_vec.push_back( tmpMat );
        }
        // Initialize a 3D matrix using the vector.
        Matrix3DXd rePart = Matrix3DXd( tmp_mat_vec );

        vector< Eigen::MatrixXd > tmp_mat_vec2;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 0+z, 0.1+z, 0.2+z, 0.3+z;
            tmp_mat_vec2.push_back( tmpMat );
        }
        // Initialize a 3D matrix using the vector.
        Matrix3DXd imPart = Matrix3DXd( tmp_mat_vec2 );

        double absRMS = Matrix3DXd::RMS_total_comp( rePart, imPart );

        double trueAns = std::pow( 1237.0/200, 0.5 );

        cout << "Absolute RMS is: " << absRMS << endl;
        cout << "Expected RMS is: " << trueAns << endl;
        if( abs( trueAns - absRMS ) < Matrix3DXd::DEF_NUM_THRESH ){
            cout << "Result: MATCH!" << endl;
        }else{
            cout << "Result: MISMATCH!" << endl;
        }
    }

}




void tests::Matrix3DXd_test_4_supp( int case_idx ){

    // Initialize test case index.
    int case_cnt = 0;

    // 0- RMS_total_comp test.
    if( case_cnt == case_idx ){

        // Create a vector of 2D matrices.
        unsigned int level_cnt = 10;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << 0+z*10, 1+z*10, 2+z*10, 3+z*10;
            tmp_mat_vec.push_back( tmpMat );
        }
        // Initialize a 3D matrix using the vector.
        Matrix3DXd myMat3D = Matrix3DXd( tmp_mat_vec );

        // Obtain subset.
        unsigned int startIdx = 2;
        unsigned int len = 3;
        Matrix3DXd myMat3Dsubset = myMat3D.segment( startIdx, len );

        bool check_sub_len = myMat3Dsubset.levels() == len;
        cout << "Is sub length correct: " << check_sub_len << endl;

        bool check_val_a_0 = myMat3Dsubset.at(0)(0,0) == 20 && 
            myMat3Dsubset.at(0)(0,1) == 21 &&
            myMat3Dsubset.at(0)(1,0) == 22 &&
            myMat3Dsubset.at(0)(1,1) == 23;
        cout << "Are values at index = 0 correct: " << check_val_a_0 << endl;

    }

}