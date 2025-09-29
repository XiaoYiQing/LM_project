#pragma once

#include <algorithm>
#include <random>
#include <vector>

using namespace std;


/* Generate a vector of random integers. */
vector<int>* randIntGen( int L_bnd, int U_bnd, unsigned int cnt );

/* Generate a vector of random integers. */
vector<int> randIntVectGen( int L_bnd, int U_bnd, unsigned int cnt );

template<typename T, typename A>
// Function to shuffle a vector
void shuffleVector(std::vector<T,A>& inVect) {
    // Obtain a random number generator
    std::random_device rd; // Obtain random number from hardware
    std::mt19937 generator(rd()); // Seed the generator

    // Shuffle the vector
    std::shuffle(inVect.begin(), inVect.end(), generator);
}

/* 
Cut-off the upper 10th order digits of the target number down to and 
including the specified order. 
*/
long long cutUpperOrder( long long tarNum, int orderPos );

/*
Turns the integer (long long) into a vector of integers where each entry is the value
of the number at a 10th power.
*/
vector<int> intToArray(long long number);

/*
Simply determines if two doubles are close enough in value with specified threshold.
*/
bool isEqEnough( double a, double b, double thresh = std::pow( 10, -9 ) );



namespace utils{

// ====================================================================== >>>>>
//      Random Number Generator
// ====================================================================== >>>>>

vector<int> rIntGen( int L_bnd, int U_bnd, unsigned int cnt );

shared_ptr<vector<double>> rDoubleGen( double L_bnd, double U_bnd, unsigned int cnt );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Indexing Arrays Functions
// ====================================================================== >>>>>

/*
Generate the index array of even indices within the specified range.
Array is in ascending order.
*/
vector< unsigned int > gen_even_idx_arr( unsigned int lower, unsigned int upper );

/*
Generate the index array of odd indices within the specified range.
Array is in ascending order.
*/
vector< unsigned int > gen_odd_idx_arr( unsigned int lower, unsigned int upper );

/*
Function generates a linearly distributed array of indices between the specified
range and filled with the specified number of indices.
Note:
    - The low and high bounds are included in the final index array.
    - There can be no repeat of indices. If this cannot be avoided given the input
        specifications, an exception is thrown.
*/
vector< unsigned int > gen_lin_idx_arr( unsigned int lower, unsigned int upper, unsigned int cnt );

/*
Functions generates an index array of ascending values between (and including) the 
specified lower and upper bounds and which ARE NOT IN the specified existing index 
array "p1".
Note:
    - "p1" must not contain repeating values.
    - "p1" must only contain integers within the specified bounds (Can be bound values).
*/
vector< unsigned int > gen_rem_idx_arr( unsigned int lower, unsigned upper, 
    vector< unsigned int >& p1 );

// ====================================================================== <<<<<


// ====================================================================== >>>>>
//      Numerical Vector Utilities
// ====================================================================== >>>>>

/*
Function generates a vector containing all matching numerical contents between
the two input numerical vectors.
*/
template <typename Z>
vector<Z> gen_match_vector( const vector<Z>& vec_A, const vector<Z>& vec_B ){
    
    // Make sure this function is only used with numerical vector.
    static_assert( std::is_arithmetic<Z>::value, "Vector contains non-numeric type" );

    vector<Z> vec_A_cp = vec_A;
    vector<Z> vec_B_cp = vec_B;

    // Sort the vectors
    sort(vec_A_cp.begin(), vec_A_cp.end());
    sort(vec_B_cp.begin(), vec_B_cp.end());

    vector<Z> matchedEntries;

    // Perform set intersection
    set_intersection(
        vec_A_cp.begin(), vec_A_cp.end(),
        vec_B_cp.begin(), vec_B_cp.end(),
        std::back_inserter( matchedEntries )
    );

    return matchedEntries;

}


/*
Sort the given numerical vector.
If "ascend" is true, sort in ascending order, else, in descending order.
*/
template <typename Z>
void sort_num_vec_inplace( vector<Z>& vec_A, bool ascend ){

    // Make sure this function is only used with numerical vector.
    static_assert( std::is_arithmetic<Z>::value, "Vector contains non-numeric type" );

    if( ascend ){
        // Sort in ascending order.
        std::sort( vec_A.begin(), vec_A.end() );
    }else{
        // Sort in descending order.
        std::sort( vec_A.begin(), vec_A.end(), std::greater<>() );
    }

}

// ====================================================================== <<<<<


}

