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
    unsigned int out_cnt = f1Data.get_out_cnt();
    unsigned int in_cnt = f1Data.get_in_cnt();

    // Determine the expected size of the Loewner Matrix.
    unsigned int row_cnt = f2Size*out_cnt;
    unsigned int col_cnt = f1Size*in_cnt;

    // Initialize the Loewner Matrix.
    Eigen::MatrixXd LM = Eigen::MatrixXd( row_cnt, col_cnt );
    // Initialize temporary matrix representing the current sub-block being computed.
    Eigen::MatrixXd currBlk = Eigen::MatrixXd( out_cnt, in_cnt );

    // Initialize current lead coordinate where the sub-block is to be inserted in the
    // Loewner Matrix.
    unsigned int lead_x = 0, lead_y = 0;

    for( unsigned int i = 0; i < f2Size; i++ ){
        for( unsigned int j = 0; j < f2Size; j++ ){

            lead_x = i*row_cnt;
            lead_y = j*col_cnt;

            

        }
    }

    return LM;
    
}
