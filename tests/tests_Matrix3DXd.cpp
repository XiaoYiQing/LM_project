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
        }catch (const std::out_of_range& e){
            cerr << e.what() << endl;
            return;
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

        Matrix3DXd my3DMat = Matrix3DXd( tmp_mat_vec );
        cout << my3DMat.at(2) << endl;
        Matrix3DXd my3DMat2 = my3DMat*3;
        cout << my3DMat2.at(2) << endl;
                
    }
    
    case_cnt++;
    // 2- Inplace multiplication operator test.
    if( case_cnt == case_idx ){

        unsigned int level_cnt = 4;
        vector< Eigen::MatrixXd > tmp_mat_vec;
        for( unsigned int z = 0; z < level_cnt; z++ ){
            Eigen::MatrixXd tmpMat = Eigen::MatrixXd( 2, 2 );
            tmpMat << z,z,z,z;
            tmp_mat_vec.push_back( tmpMat );
        }

        Matrix3DXd my3DMat = Matrix3DXd( tmp_mat_vec );
        cout << my3DMat.at(2) << endl;
        my3DMat *= 3;
        cout << my3DMat.at(2) << endl;
                
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

        try{
            Matrix3DXd my_3D_mat = Matrix3DXd( tmp_mat_vec );
        }catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
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

        // Add another matrix to the vector at its end.
        Eigen::MatrixXd addMat = Eigen::MatrixXd( 2, 2 );
        addMat << 91, 92, 93, 94;
        try{
            my_3D_mat.push_back( addMat );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        }
        cout << my_3D_mat.at( my_3D_mat.levels() - 1 ) << endl;

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
        cout << my_3D_mat.at( tarIdx ) << endl;
        cout << my_3D_mat.at( tarIdx - 1 ) << endl;


        // Insert another matrix to the vector at the target invalid index.
        Eigen::MatrixXd addMat2 = Eigen::MatrixXd( 2, 2 );
        addMat2 << 51, 52, 53, 54;
        tarIdx = 100;
        try{
            my_3D_mat.insert( tarIdx, addMat2 );
        } catch (const std::invalid_argument& e){
            cerr << e.what() << endl;
            return;
        } catch( const std::out_of_range& e ){
            cerr << e.what() << endl;
            return;
        }
        cout << my_3D_mat.at( tarIdx ) << endl;
        cout << my_3D_mat.at( tarIdx - 1 ) << endl;

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

}
