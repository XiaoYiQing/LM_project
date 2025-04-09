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

