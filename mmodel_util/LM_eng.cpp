#include "LM_eng.h"





shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_LM( const fData& f1Data, const fData& f2Data ){

    // Obtain the size of the two partitions.
    unsigned int f1Size = f1Data.get_f_cnt();
    unsigned int f2Size = f2Data.get_f_cnt();
    if( f1Size == 0 || f2Size == 0 ){
        throw std::invalid_argument( "Empty frequency data inputs are not allowed for constructing the Loewner Matrix." );
    }
    
    // Obtain the number of outputs and inputs.
    unsigned int out_cnt = f1Data.get_out_cnt();
    unsigned int in_cnt = f1Data.get_in_cnt();

    // Determine the expected size of the Loewner Matrix.
    unsigned int row_cnt = f2Size*out_cnt;
    unsigned int col_cnt = f1Size*in_cnt;

    // Initialize the Loewner Matrix.
    shared_ptr<Eigen::MatrixXcd> LM = std::make_shared<Eigen::MatrixXcd>( row_cnt, col_cnt );

    // Initialize temporary matrix representing the current sub-block being computed.
    Eigen::MatrixXcd LM_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Initialize current freq. data from each partition.
    Eigen::MatrixXcd f1_D_j = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd f2_D_i = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Initialize current frequencies from each partition.
    complex<double> f1_j = complex<double>(0,0); 
    complex<double> f2_i = complex<double>(0,0);
    

    // Initialize current lead coordinate where the sub-block is to be inserted in the
    // Loewner Matrix.
    unsigned int lead_x = 0, lead_y = 0;

    // Build the Loewner Matrix block by block.
    for( unsigned int i = 0; i < f2Size; i++ ){

        // Obtain the current f2 freq. and data.
        f2_i = f2Data.get_cplx_f_at( i );
        f2_D_i = f2Data.get_cplxData_at_f( i );

        for( unsigned int j = 0; j < f1Size; j++ ){

            // Obtain the current f1 freq. and data.
            f1_j = f1Data.get_cplx_f_at( j );
            f1_D_j = f1Data.get_cplxData_at_f( j );

            // Obtain the matrix coordinate of the upper-left leading point of the current
            // LM sub-block being computed.
            lead_x = i*out_cnt;
            lead_y = j*in_cnt;

            // Compute the current Loewner Matrix block.
            LM_ij = ( f2_D_i - f1_D_j )/( f2_i - f1_j );

            // Insert the current calculated block into its part in the full matrix.
            LM->block( lead_x, lead_y, out_cnt, in_cnt) = LM_ij;

        }

    }


    return LM;

}



shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_SLM( const fData& f1Data, const fData& f2Data ){

    // Obtain the size of the two partitions.
    unsigned int f1Size = f1Data.get_f_cnt();
    unsigned int f2Size = f2Data.get_f_cnt();
    if( f1Size == 0 || f2Size == 0 ){
        throw std::invalid_argument( "Empty frequency data inputs are not allowed for constructing the shifted Loewner Matrix." );
    }

    // Obtain the number of outputs and inputs.
    unsigned int out_cnt = f1Data.get_out_cnt();
    unsigned int in_cnt = f1Data.get_in_cnt();

    // Determine the expected size of the shifted Loewner Matrix.
    unsigned int row_cnt = f2Size*out_cnt;
    unsigned int col_cnt = f1Size*in_cnt;

    // Initialize the shifted Loewner Matrix.
    shared_ptr<Eigen::MatrixXcd> sLM = std::make_shared<Eigen::MatrixXcd>( row_cnt, col_cnt );

    // Initialize temporary matrix representing the current sub-block being computed.
    Eigen::MatrixXcd sLM_ij = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Initialize current freq. data from each partition.
    Eigen::MatrixXcd f1_D_j = Eigen::MatrixXcd( out_cnt, in_cnt );
    Eigen::MatrixXcd f2_D_i = Eigen::MatrixXcd( out_cnt, in_cnt );
    // Initialize current frequencies from each partition.
    complex<double> f1_j = complex<double>(0,0); 
    complex<double> f2_i = complex<double>(0,0);

    // Initialize current lead coordinate where the sub-block is to be inserted in the
    // shifted Loewner Matrix.
    unsigned int lead_x = 0, lead_y = 0;

    // Build the shifted Loewner Matrix block by block.
    for( unsigned int i = 0; i < f2Size; i++ ){

        // Obtain the current f2 freq. and data.
        f2_i = f2Data.get_cplx_f_at( i );
        f2_D_i = f2Data.get_cplxData_at_f( i );

        for( unsigned int j = 0; j < f1Size; j++ ){

            // Obtain the current f1 freq. and data.
            f1_j = f1Data.get_cplx_f_at( j );
            f1_D_j = f1Data.get_cplxData_at_f( j );

            // Obtain the matrix coordinate of the upper-left leading point of the current
            // sLM sub-block being computed.
            lead_x = i*out_cnt;
            lead_y = j*in_cnt;

            // Compute the current shifted Loewner Matrix block.
            sLM_ij = ( f2_i*f2_D_i - f1_j*f1_D_j )/( f2_i - f1_j );

            // Insert the current calculated block into its part in the full matrix.
            sLM->block( lead_x, lead_y, out_cnt, in_cnt) = sLM_ij;

        }

    }

    return sLM;

}

shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_W( const fData& f1Data ){

    // Obtain the size of the two partitions.
    unsigned int f1Size = f1Data.get_f_cnt();
    if( f1Size == 0 ){
        throw std::invalid_argument( "Empty frequency data input is not allowed for constructing the W matrix vector." );
    }

    // Obtain the number of outputs and inputs.
    unsigned int out_cnt = f1Data.get_out_cnt();
    unsigned int in_cnt = f1Data.get_in_cnt();

    // Determine the expected size of the shifted Loewner Matrix.
    unsigned int col_cnt = f1Size*in_cnt;

    // Initialize the W matrix vector.
    shared_ptr<Eigen::MatrixXcd> W = std::make_shared<Eigen::MatrixXcd>( out_cnt, col_cnt );

    // Initialize current freq. data from each partition.
    Eigen::MatrixXcd f1_D_j = Eigen::MatrixXcd( out_cnt, in_cnt );

    // Initialize current lead coordinate where the sub-block is to be inserted in the
    // shifted Loewner Matrix.
    unsigned int lead_x = 0, lead_y = 0;

    for( unsigned int j = 0; j < f1Size; j++ ){

        // Obtain the current f1 freq. and data.
        f1_D_j = f1Data.get_cplxData_at_f( j );

        // Obtain the matrix coordinate of the upper-left leading point of the current
        // sLM sub-block being computed.
        lead_y = j*in_cnt;

        // Insert the current calculated block into its part in the full matrix.
        W->block( lead_x, lead_y, out_cnt, in_cnt) = f1_D_j;

    }

    return W;

}


shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_F( const fData& f2Data ){

    // Obtain the size of the two partitions.
    unsigned int f2Size = f2Data.get_f_cnt();
    if( f2Size == 0 ){
        throw std::invalid_argument( "Empty frequency data input is not allowed for constructing the F matrix vector." );
    }

    // Obtain the number of outputs and inputs.
    unsigned int out_cnt = f2Data.get_out_cnt();
    unsigned int in_cnt = f2Data.get_in_cnt();

    // Determine the expected size of the F matrix vector.
    unsigned int row_cnt = f2Size*out_cnt;

    // Initialize the W matrix vector.
    shared_ptr<Eigen::MatrixXcd> F = std::make_shared<Eigen::MatrixXcd>( row_cnt, in_cnt );

    // Initialize current freq. data from each partition.
    Eigen::MatrixXcd f2_D_i = Eigen::MatrixXcd( out_cnt, in_cnt );

    // Initialize current lead coordinate where the sub-block is to be inserted in the
    // shifted Loewner Matrix.
    unsigned int lead_x = 0, lead_y = 0;

    for( unsigned int i = 0; i < f2Size; i++ ){

        // Obtain the current f1 freq. and data.
        f2_D_i = f2Data.get_cplxData_at_f( i );

        // Obtain the matrix coordinate of the upper-left leading point of the current
        // sLM sub-block being computed.
        lead_x = i*out_cnt;

        // Insert the current calculated block into its part in the full matrix.
        F->block( lead_x, lead_y, out_cnt, in_cnt) = f2_D_i;

    }

    return F;

}


shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_LM_pencil( complex<double> ref_f, const Eigen::MatrixXcd& LM, const Eigen::MatrixXcd& SLM ){


    unsigned int row_cnt = LM.rows();
    unsigned int col_cnt = LM.cols();
    if( row_cnt != SLM.rows() || col_cnt != SLM.cols() ){
        throw std::invalid_argument( "The LM and the SLM must shared the same dimensions." );
    }
    
    shared_ptr<Eigen::MatrixXcd> LM_pen = std::make_shared<Eigen::MatrixXcd>( row_cnt, col_cnt );

    *LM_pen = LM - ref_f*SLM;

    return LM_pen;

}


shared_ptr<Eigen::MatrixXcd> LM_UTIL::build_reT_mat( bool has_DC_pt, unsigned int sub_mat_size, unsigned int sub_blk_cnt ){


    // Obtain the total number of rows and columns of the final transformation matrix.
    unsigned int row_cnt = 2*sub_blk_cnt*sub_mat_size;
    unsigned int col_cnt = 2*sub_blk_cnt*sub_mat_size;

    if( has_DC_pt ){
        row_cnt += sub_mat_size;
        col_cnt += sub_mat_size;        
    }

    // Complex 1.
    complex<double> im1(0,1);
    // Initialize the standard identity matrix.
    Eigen::MatrixXd I_mat = Eigen::MatrixXd::Identity(sub_mat_size, sub_mat_size);
    // Initialize the unit sub-block transformation matrix.
    Eigen::MatrixXcd T_unit = Eigen::MatrixXcd( 2*sub_mat_size, 2*sub_mat_size );
    T_unit.block( 0, 0, sub_mat_size, sub_mat_size ) = I_mat;
    T_unit.block( 0, sub_mat_size, sub_mat_size, sub_mat_size ) = -im1*I_mat;
    T_unit.block( sub_mat_size, 0, sub_mat_size, sub_mat_size ) = I_mat;
    T_unit.block( sub_mat_size, sub_mat_size, sub_mat_size, sub_mat_size ) = im1*I_mat;

    // Initialize the transformation matrix.
    shared_ptr<Eigen::MatrixXcd> T_mat = std::make_shared<Eigen::MatrixXcd>( row_cnt, col_cnt );
    T_mat->setZero();
    
    // Initialize current sub-block's lead coordinate.
    unsigned int lead_x = 0;
    if( has_DC_pt ){
        T_mat->block( lead_x, lead_x, sub_mat_size, sub_mat_size ) = I_mat;
        lead_x += sub_mat_size;
    }

    for( unsigned int z = 0; z < sub_blk_cnt; z++ ){

        // Insert the current sub-block.
        T_mat->block( lead_x, lead_x, 2*sub_mat_size, 2*sub_mat_size ) = T_unit;

        // Increment leading point on the matrix diagonal.
        lead_x += 2*sub_mat_size;

    }

    return T_mat;

}

