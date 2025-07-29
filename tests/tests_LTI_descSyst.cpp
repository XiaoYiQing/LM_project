#include "tests_LTI_descSyst.h"




void tests::LTI_descSyst_test_1( unsigned int case_idx ){


    // Initialize test case index.
    int case_cnt = 0;

    // 0- Test base functionalities such as initialization, insert, get.
    if( case_cnt == case_idx ){

        // Number of inputs.
        unsigned int m = 7;
        
        // Number of outputs.
        unsigned int p = 5;

        // System order
        unsigned int n = 40;
        
        shared_ptr< Eigen::MatrixXd >E_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Random( n, n ) );
        shared_ptr< Eigen::MatrixXd >A_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Random( n, n ) );
        shared_ptr< Eigen::MatrixXd >B_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Random( n, p ) );
        shared_ptr< Eigen::MatrixXd >C_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Random( m, n ) );
        shared_ptr< Eigen::MatrixXd >D_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Random( m, p ) );


        LTI_descSyst mySyst = LTI_descSyst();
        mySyst.set_E( E_ptr );
        if( mySyst.is_consistent() ){
            cout << "Consistency should be false!" << endl;
            return;
        }
        mySyst.set_A( A_ptr );
        if( mySyst.is_consistent() ){
            cout << "Consistency should be false!" << endl;
            return;
        }
        mySyst.set_B( B_ptr );
        if( mySyst.is_consistent() ){
            cout << "Consistency should be false!" << endl;
            return;
        }
        mySyst.set_C( C_ptr );
        if( mySyst.is_consistent() ){
            cout << "Consistency should be false!" << endl;
            return;
        }
        mySyst.set_D( D_ptr );
        if( !mySyst.is_consistent() ){
            cout << "Consistency should be true!" << endl;
            return;
        }

        if( mySyst.get_input_cnt() != m ){
            cout << "Input count mismatch." << endl;
            return;
        }
        if( mySyst.get_output_cnt() != p ){
            cout << "Output count mismatch." << endl;
            return;
        }
        if( mySyst.get_order() != n ){
            cout << "Order mismatch." << endl;
            return;
        }

        if( mySyst.get_E() != *E_ptr ){
            cout << "Matrix E compare failed." << endl;
            return;
        }
        if( mySyst.get_A() != *A_ptr ){
            cout << "Matrix A compare failed." << endl;
            return;
        }
        if( mySyst.get_B() != *B_ptr ){
            cout << "Matrix B compare failed." << endl;
            return;
        }
        if( mySyst.get_C() != *C_ptr ){
            cout << "Matrix C compare failed." << endl;
            return;
        }
        if( mySyst.get_D() != *D_ptr ){
            cout << "Matrix D compare failed." << endl;
            return;
        }

        cout << "All tests passed!" << endl;

    }

    case_cnt++;
    // 1- Test stability checking.
    if( case_cnt == case_idx ){

        // Number of inputs.
        unsigned int m = 7;
        
        // Number of outputs.
        unsigned int p = 5;

        // System order
        unsigned int n = 100;
        
        shared_ptr< Eigen::MatrixXd >E_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Identity( n, n ) );
        shared_ptr< Eigen::MatrixXd >A_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Identity( n, n ) );
        shared_ptr< Eigen::MatrixXd >B_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Random( n, p ) );
        shared_ptr< Eigen::MatrixXd >C_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Random( m, n ) );
        shared_ptr< Eigen::MatrixXd >D_ptr = make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Random( m, p ) );

        LTI_descSyst mySyst = LTI_descSyst();
        mySyst.set_E( E_ptr );
        mySyst.set_A( A_ptr );
        mySyst.set_B( B_ptr );
        mySyst.set_C( C_ptr );
        mySyst.set_D( D_ptr );
        bool is_stab = mySyst.is_stable();

        if( is_stab ){
            cout << "Failed test: found stable but case was not stable." << endl;
            return;
        }else{
            cout << "Stability check 1: passed" << endl;
        }

        *A_ptr = -1*Eigen::MatrixXd::Identity( n, n );
        mySyst.set_A( A_ptr );
        is_stab = mySyst.is_stable();
        if( !is_stab ){
            cout << "Failed test: found not stable but case was stable." << endl;
            return;
        }else{
            cout << "Stability check 2: passed" << endl;
        }

        *A_ptr = -1*Eigen::MatrixXd::Identity( n-1, n-1 );
        mySyst.set_A( A_ptr );
        is_stab = mySyst.is_stable();
        if( is_stab ){
            cout << "Failed test: inconsistent system found stable." << endl;
            return;
        }else{
            cout << "Stability check 3: passed" << endl;
        }

    }

    case_cnt++;
    // 2- Transfer function eval test.
    if( case_cnt == case_idx ){

        // Number of inputs.
        unsigned int m = 2;
        
        // Number of outputs.
        unsigned int p = 2;

        // System order
        unsigned int n = 4;
        
        shared_ptr< Eigen::MatrixXd >E_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Identity( n, n ) );
        shared_ptr< Eigen::MatrixXd >A_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Identity( n, n ) );
        shared_ptr< Eigen::MatrixXd >B_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Ones( n, p ) );
        shared_ptr< Eigen::MatrixXd >C_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Ones( m, n ) );
        shared_ptr< Eigen::MatrixXd >D_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Zero( m, p ) );

        LTI_descSyst mySyst = LTI_descSyst();
        mySyst.set_E( E_ptr );
        mySyst.set_A( A_ptr );
        mySyst.set_B( B_ptr );
        mySyst.set_C( C_ptr );
        mySyst.set_D( D_ptr );

        double test_val = 0;
        Eigen::MatrixXcd eval_mat = mySyst.tf_eval( test_val );
        Eigen::MatrixXcd expect_mat = -4*Eigen::MatrixXcd::Ones( m, p );

        if( eval_mat == expect_mat ){
            cout << "Test 1 passed: single f evaluation." << endl;
        }else{
            cout << "Test 1 failed: single f evaluation." << endl;
        }


        vector< complex<double> > test_arr = { complex<double>(0,0), complex<double>(0.5,0),
            complex<double>(0,0.5) };
        Matrix3DXcd eval_arr = mySyst.tf_eval( test_arr );

        Matrix3DXcd expect_arr = Matrix3DXcd( m, p, test_arr.size() );

        Eigen::MatrixXcd tmp = (Eigen::MatrixXcd(2, 2) << 
        std::complex<double>(-4,0), std::complex<double>(-4,0),
        std::complex<double>(-4,0), std::complex<double>(-4,0)).finished();
        expect_arr.set( 0, tmp );
        tmp = (Eigen::MatrixXcd(2, 2) << 
        std::complex<double>(-8,0), std::complex<double>(-8,0),
        std::complex<double>(-8,0), std::complex<double>(-8,0)).finished();
        expect_arr.set( 1, tmp );
        tmp = (Eigen::MatrixXcd(2, 2) << 
        std::complex<double>(-3.2,-1.6), std::complex<double>(-3.2,-1.6),
        std::complex<double>(-3.2,-1.6), std::complex<double>(-3.2,-1.6)).finished();
        expect_arr.set( 2, tmp );

        bool test_pass = true;
        for( unsigned int z = 0; z < test_arr.size(); z++ ){
            test_pass = test_pass && ( eval_arr.at(z) == expect_arr.at(z) );
        }
        if( test_pass ){
            cout << "Test 2 passed: multi f evaluation." << endl;
        }else{
            cout << "Test 2 failed: multi f evaluation." << endl;
        }

    }



    case_cnt++;
    // 3- Saved compute poles test.
    if( case_cnt == case_idx ){

        // Number of inputs.
        unsigned int m = 2;
        // Number of outputs.
        unsigned int p = 2;
        // System order
        unsigned int n = 4;
        
        shared_ptr< Eigen::MatrixXd >E_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Identity( n, n ) );
        shared_ptr< Eigen::MatrixXd >A_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Identity( n, n ) );
        shared_ptr< Eigen::MatrixXd >B_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Ones( n, p ) );
        shared_ptr< Eigen::MatrixXd >C_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Ones( m, n ) );
        shared_ptr< Eigen::MatrixXd >D_ptr = 
            make_shared< Eigen::MatrixXd >( Eigen::MatrixXd::Zero( m, p ) );

        LTI_descSyst mySyst = LTI_descSyst();
        mySyst.set_E( E_ptr );
        mySyst.set_A( A_ptr );
        mySyst.set_B( B_ptr );
        mySyst.set_C( C_ptr );
        mySyst.set_D( D_ptr );

        mySyst.get_poles();

    }

}
