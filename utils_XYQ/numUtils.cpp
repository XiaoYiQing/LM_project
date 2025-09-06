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
//      Random Number Generator
// ====================================================================== >>>>>

shared_ptr< vector<int> > utils::rIntGen( int L_bnd, int U_bnd, unsigned int cnt ){

    if( L_bnd > U_bnd ){
        throw std::invalid_argument( "rIntGen lower bound must not be larger than the upper bound." );
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
    shared_ptr<vector<int>> tmpVec = make_shared< vector<int> >();
    tmpVec->reserve( cnt );
    for (unsigned int i = 0; i < cnt; i++) {
        tmpVec->emplace_back( distrib(gen) );
    }
    
    return tmpVec;

}

shared_ptr<vector<double>> utils::rDoubleGen( double L_bnd, double U_bnd, unsigned int cnt ){

    if( L_bnd > U_bnd ){
        throw std::invalid_argument( "rDoubleGen lower bound must not be larger than the upper bound." );
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
    std::uniform_real_distribution<double> distrib( L_bnd, U_bnd );

    // // Create an integer array to store all random integers.
    shared_ptr<vector<double>> tmpVec = make_shared< vector<double> >();
    tmpVec->reserve( cnt );
    for (unsigned int i = 0; i < cnt; i++) {
        tmpVec->emplace_back( distrib(gen) );
    }
    
    return tmpVec;

}

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Indexing Arrays Functions
// ====================================================================== >>>>>

vector< unsigned int > utils::gen_even_idx_arr( unsigned int lower, unsigned int upper ){    

    if( upper <= lower ){
        throw std::invalid_argument( "The upper bound must be larger than the lower bound." );
    }

    // Initialize return vector.
    vector< unsigned int > resVec;
    resVec.reserve( ( upper - lower )/2 + 1 );

    // Initialize the starting index.
    unsigned int start_idx = lower;
    if( remainder( (double) lower, 2.0 ) == 1 ){
        start_idx += 1;
    }

    // Fill even entries into the index array.
    for( unsigned int z = start_idx; z <= upper; z += 2 ){
        resVec.push_back( z );
    }
    resVec.shrink_to_fit();

    return resVec;

}

vector< unsigned int > utils::gen_odd_idx_arr( unsigned int lower, unsigned int upper ){

    if( upper <= lower ){
        throw std::invalid_argument( "The upper bound must be larger than the lower bound." );
    }

    // Initialize return vector.
    vector< unsigned int > resVec;
    resVec.reserve( ( upper - lower )/2 + 1 );

    // Initialize the starting index.
    unsigned int start_idx = lower;
    if( remainder( (double) lower, 2.0 ) == 0 ){
        start_idx += 1;
    }

    // Fill even entries into the index array.
    for( unsigned int z = start_idx; z <= upper; z += 2 ){
        resVec.push_back( z );
    }
    resVec.shrink_to_fit();

    return resVec;

}

vector< unsigned int > utils::gen_lin_idx_arr( unsigned int lower, unsigned int upper, unsigned int cnt ){

    if( cnt < 2 ){
        throw std::invalid_argument( "The length of the linear index array cannot be less than 2." );
    }
    if( upper <= lower ){
        throw std::invalid_argument( "The upper bound must be larger than the lower bound." );
    }
    if( upper - lower + 1 < cnt ){
        throw std::invalid_argument( "The number of indices required must not exceed the number of distinct integers between the lower and upper bounds." );
    }
    if( cnt == 2 ){
        vector< unsigned int > resVec = {lower, upper};
        return resVec;
    }

    // Compute the average distance between each index (if allowed to be floating numbers).
    double avg_diff = ( (double) ( upper - lower ) )/( cnt - 1 );
    // Initialize return vector.
    vector< unsigned int > resVec;
    resVec.reserve( cnt );

    // Add the lower bound.
    resVec.push_back( lower );
    // Add all indices between the bounds.
    for( unsigned int z = 0; z < cnt - 2; z++ ){
        resVec.push_back( (unsigned int) std::ceil( avg_diff*(z+1) + lower ) );
    }
    // Add the upper bound.
    resVec.push_back( upper );

    return resVec;

}


vector< unsigned int > utils::gen_rem_idx_arr( unsigned int lower, unsigned upper, vector< unsigned int >& p1 ){

    vector< unsigned int > p1_copy = p1;
    // Sort the vector in ascending order
    std::sort(p1_copy.begin(), p1_copy.end());

    // Index outside range check.
    if( *( p1_copy.end() - 1 ) > upper ){
        throw std::out_of_range( "The given index array \"p1\" contains at least one entry outside the specified range." );
    }
    // Add an outside range value at the end of the reference vector to serve as
    // an impossible value to reach.
    p1_copy.push_back( upper + 1 );

    // Initialize return index vector.
    vector< unsigned int > retVec;
    retVec.reserve( upper - lower + 1 );

    // Initialize variables 
    unsigned int retVec_idx = 0;
    unsigned int p1_idx = 0;
    
    // Cycle through all possible indices within the specified range and insert
    // into the result vector if conditions are met.
    for( unsigned int currIdx = lower; currIdx <= upper; currIdx++ ){

        if( currIdx == p1_copy.at( p1_idx ) ){
            p1_idx++;
        }else{
            retVec.push_back( currIdx );
        }

    }

    return retVec;
    
}


// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Numerical Vector Utilities
// ====================================================================== >>>>>

template <typename T>
vector<T> utils::gen_match_vector( vector<T> vec_A, vector<T> vec_B ){
    
    static_assert(std::is_arithmetic<T>::value, "Template parameter must be a numerical type");

    vector<T> matchedEntries;

    // Perform set intersection
    std::set_intersection(
        vec_A.begin(), vec_A.end(),
        vec_B.begin(), vec_B.end(),
        std::back_inserter( matchedEntries )
    );

    return matchedEntries;
    
}

// ====================================================================== <<<<<

