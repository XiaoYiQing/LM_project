#include "tests_Eigen.h"


using namespace std;

/*
This is simply a function where I test random things to learn about how to use 
the eigen library.
*/
void tests::eigen_test1( int testCase ){


    int case_cnt = 0;

// ---------------------------------------------------------------------- >>>>>
//      Simple 2x2 Matrix Definition
// ---------------------------------------------------------------------- >>>>>
    if( case_cnt == testCase ){       
        Eigen::Matrix2d matA;
        matA << 1, 2, 3, 4;
        Eigen::Matrix2d matA_inv = matA.inverse();
        cout << "Matrix A:\n" << matA << endl;
        cout << "Matrix A inverse:\n" << matA_inv << endl;
    }

    case_cnt++;
// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      How to Increase Vector Size
// ---------------------------------------------------------------------- >>>>>

    if( case_cnt == testCase ){
        // How to add entries to a vector.
        Eigen::VectorXd vecA(2);
        vecA << 1, 2;
        cout << "Vector A:\n" << vecA << endl;
        // Resize the vector to the new size
        vecA.conservativeResize(5);
        vecA(2) = 3;    vecA(3) = 4;    vecA(4) = 55;
        cout << "Vector A after size increase:\n" << vecA << endl;
    }

    case_cnt++;
// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      How to Increase 2D Matrix Size
// ---------------------------------------------------------------------- >>>>>

    if( case_cnt == testCase ){
        // A 2D matrix of variable size.
        int row_cnt = 4; // Rows
        int col_cnt = 5; // Columns
        Eigen::MatrixXd matX( row_cnt, col_cnt );
        for( int i = 0; i < row_cnt; i++ ){
            for( int j = 0; j < col_cnt; j++ ){
                matX( i, j ) = i+j;
            }
        }
        cout << "Matrix X:\n" << matX << endl;

        // Increase the size to 5x7
        int row_cnt2 = 5; // Rows
        int col_cnt2 = 7; // Columns
        // Increase size. NOTE: if specified size is less than original, you simply
        // lose the data outside the new size range.
        matX.conservativeResize(5, 7);
        for( int i = 0; i < row_cnt2; i++ ){
            for( int j = col_cnt; j < col_cnt2; j++ ){
                matX( i, j ) = i+j;
            }
        }
        for( int i = row_cnt; i < row_cnt2; i++ ){
            for( int j = 0; j < col_cnt; j++ ){
                matX( i, j ) = i+j;
            }
        }
        cout << "Matrix X:\n" << matX << endl;
    }

    case_cnt++;
// ---------------------------------------------------------------------- <<<<<


// ---------------------------------------------------------------------- >>>>>
//      How to Create a 3D Matrix
// ---------------------------------------------------------------------- >>>>>
    
    if( case_cnt == testCase ){
        // Define the number of 2D matrix in our 3D matrix.
        int mat_cnt = 3;
        int row_cnt = 4;
        int col_cnt = 5;

        /*
        There is no direct 3D matrix representation in Eigen.
        You have to create a collection of 2D matrices.
        */
        vector<Eigen::MatrixXd> threeDMatrix( mat_cnt, Eigen::MatrixXd( row_cnt, col_cnt ) );

        // Initialize the 3D matrix with some values
        for (int z = 0; z < mat_cnt; z++) {
            for (int i = 0; i < row_cnt; i++) {
                for (int j = 0; j < col_cnt; j++) {
                    threeDMatrix[z](i, j) = z * 10 + i + j; // Example initialization
                }
            }
        }

        // Output the 3D matrix values
        for (int k = 0; k < mat_cnt; k++) {
            cout << "Slice " << k << ":\n" << threeDMatrix[k] << "\n";
        }

        int lol = 0;
    }

    case_cnt++;
// ---------------------------------------------------------------------- <<<<<

}


