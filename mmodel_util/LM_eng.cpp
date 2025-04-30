#include "LM_eng.h"



Eigen::MatrixXd LM_UTIL::build_LM( fData f1Data, fData f2Data ){

    // Obtain the size of the two partitions.
    unsigned int f1Size = f1Data.get_f_cnt();
    unsigned int f2Size = f2Data.get_f_cnt();
    if( f1Size == 0 || f2Size == 0 ){
        throw std::invalid_argument( "Empty frequency data inputs are not allowed for constructing the Loewner Matrix." );
    }

    Eigen::MatrixXd refMat = f1Data.get_reData_at_f(0);
    
    // Obtain the number of outputs and inputs.
    int out_cnt = f1Data.get_out_cnt();
    int in_cnt = f1Data.get_in_cnt();

    // Determine the expected size of the Loewner Matrix.
    int row_cnt = f2Size*out_cnt;
    int col_cnt = f1Size*in_cnt;

    // // Initialize the Loewner Matrix.
    // Eigen::MatrixXcd LM = Eigen::MatrixXcd( row_cnt, col_cnt );
    // // Initialize temporary matrix representing the current sub-block being computed.
    // Eigen::MatrixXcd LM_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    // // Initialize current freq. data from each partition.
    // Eigen::MatrixXcd f1_D_j = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Eigen::MatrixXcd f2_D_i = Eigen::MatrixXcd( out_cnt, in_cnt );
    // // Initialize current frequencies from each partition.
    // complex<double> f1_j = 0, f2_i = 0;
    

    // // Initialize current lead coordinate where the sub-block is to be inserted in the
    // // Loewner Matrix.
    // unsigned int lead_x = 0, lead_y = 0;

    // // Build the Loewner Matrix block by block.
    // for( unsigned int i = 0; i < f2Size; i++ ){

    //     // // Obtain the current f2 freq. and data.
    //     // f2_i = f2Data.get_cplx_f_at( i );
    //     // f2_D_i = f2Data.get_cplxData_at_f( i );

    //     for( unsigned int j = 0; j < f1Size; j++ ){

    //         // // Obtain the current f1 freq. and data.
    //         // f1_j = f1Data.get_cplx_f_at( j );
    //         // f1_D_j = f1Data.get_cplxData_at_f( j );

    //         // // Obtain the matrix coordinate of the upper-left leading point of the current
    //         // // LM sub-block being computed.
    //         // lead_x = i*row_cnt;
    //         // lead_y = j*col_cnt;

    //         // // Compute the current Loewner Matrix block.
    //         // LM_ij = ( f2_D_i - f1_D_j )/( f2_i - f1_j );

    //         // // Insert the current calculated block into its part in the full matrix.
    //         // LM.block( lead_x, lead_y, out_cnt, in_cnt) = LM_ij;

    //     }

    // }

    // return LM;

}
