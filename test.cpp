
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <random>

#include "rref.hpp"


int main( int argc, char **argv, char **envp ) {

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(-100,100);

	
//  1    2   -1   -4
//  2    3   -1   -11
// -2    0   -3    22

    int M = 3 ;
    int N = 5 ;
    double data[] = { -3,1,2,   6,-2,-4,  -1,2,5,  1,3,8,  -7,-1,-4 } ;
    //double data[] = { 2,-3,1,1,4,1,7,-5,4,-7,-6,-5,2,3,2 } ;
    // double *data = (double*)std::malloc( sizeof(double) * M * N ) ;
    // for( int i=0 ; i<M*N ; i++ ) {
    //     data[i] = distribution(generator) ;
    // }

    matprint( data, M, N ) ;
    int freeVars = nullspace( data, M, N ) ;
    matprint( data, N, freeVars ) ;

    return 0 ;
}


