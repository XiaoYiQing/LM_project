#include "numUtils.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

vector<int>* randIntGen( int L_bnd, int U_bnd, unsigned int cnt ){

    try{
        if( L_bnd > U_bnd ){
            string errMsg = "randIntGen lower bound must not be larger than the upper bound."; 
            throw errMsg;
        }
    }catch( string errMsg ){
        cout << errMsg << endl;
        vector<int>* tmpVec = new vector<int>();
        return tmpVec;
    }

    /*
    A random number seed from hardware which ensures a varying seed based 
    on external entropy.
    */
    random_device rd;
    /*
    The generator using the seed. 
    mt19937 is a popular choice due to its balance of speed and randomness quality
    */
    mt19937 gen(rd());


    // Define a uniform distribution for even number occurence chance over the range.
    uniform_int_distribution<> distrib( L_bnd, U_bnd );

    // Create an integer array to store all random integers.
    vector<int>* tmpVec = new vector<int>( cnt );
    for (unsigned int i = 0; i < cnt; i++) {
        tmpVec->at(i) = distrib(gen);
    }

    return tmpVec;
}

vector<int> randIntVectGen( int L_bnd, int U_bnd, unsigned int cnt ){

    try{
        if( L_bnd > U_bnd ){
            string errMsg = "randIntGen lower bound must not be larger than the upper bound."; 
            throw errMsg;
        }
    }catch( string errMsg ){
        cout << errMsg << endl;
        vector<int> tmpVec;
        return tmpVec;
    }

    /*
    A random number seed from hardware which ensures a varying seed based 
    on external entropy.
    */
    random_device rd;
    /*
    The generator using the seed. 
    mt19937 is a popular choice due to its balance of speed and randomness quality
    */
    mt19937 gen(rd());


    // Define a uniform distribution for even number occurence chance over the range.
    uniform_int_distribution<> distrib( L_bnd, U_bnd );

    // Create an integer array to store all random integers.
    vector<int> tmpVec;
    for (unsigned int i = 0; i < cnt; i++) {
        tmpVec.push_back( distrib( gen ) );
    }

    return tmpVec;

}

// template<typename T, typename A>
// // Function to shuffle a vector
// void shuffleVector(std::vector<T,A>& inVect) {
//     // Obtain a random number generator
//     std::random_device rd; // Obtain random number from hardware
//     std::mt19937 generator(rd()); // Seed the generator

//     // Shuffle the vector
//     std::shuffle(inVect.begin(), inVect.end(), generator);
// }

/*
    Retains portion of the target integer "tarNum" up to decimal 
    position "orderPos".
    For example, 
        cutUpperOrder( 987654321, 5 ) => 54321
*/
long long cutUpperOrder( long long tarNum, int orderPos ){

    int rem_i = 0;
    long long cutNum = 0;
    long long tmp = 0;

    for( int i = 1; i <= orderPos; i++ ){
        rem_i = tarNum % 10;
        tmp = ( long long ) pow( 10, i - 1 );
        cutNum = cutNum + rem_i*tmp;

        tarNum = tarNum/10;
    }

    return cutNum;

}



vector<int> intToArray(long long number) {

    // TODO: make sure the vector is created on the heap.
    vector<int> digits;
    if (number == 0) {
        digits.push_back(0);
        return digits;
    }

    while (number > 0) {
        int digit = number % 10; // Extract the last digit
        digits.push_back(digit); 
        number /= 10; // Remove the last digit
    }
    // The digits are in reverse order, so reverse them
    reverse(digits.begin(), digits.end()); 
    return digits;

}



bool isEqEnough( double a, double b, double thresh ){
    return std::abs( a - b ) < std::abs( thresh );
}




// ====================================================================== >>>>>
//      Indexing Arrays Functions
// ====================================================================== >>>>>

vector< unsigned int > utils::gen_lin_idx_arr( unsigned int lower, unsigned int upper, unsigned int cnt ){

    if( cnt < 2 ){
        throw std::invalid_argument( "The length of the linear index array cannot be less than 2." );
    }
    if( upper <= lower ){
        throw std::invalid_argument( "The upper bound must be larger than the lower bound." );
    }
    if( cnt == 2 ){
        vector< unsigned int > resVec = {lower, upper};
        return resVec;
    }
    

    // Compute the average distance between 
    double avg_diff = ( (double) ( upper - lower ) )/( cnt - 1 );

    vector< unsigned int > resVec;

    resVec.push_back( lower );
    for( unsigned int z = 0; z < cnt - 2; z++ ){

        double tmp = avg_diff*(z+1) + lower;
        resVec.push_back( std::ceil( tmp ) );

    }
    resVec.push_back( upper );

    return resVec;

}


// ====================================================================== <<<<<