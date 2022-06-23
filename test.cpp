
#include <iostream>
#include "rref.hpp"

void matprint( double *data, int M, int N ) ;


int main( int argc, char **argv, char **envp ) {

//  1    2   -1   -4
//  2    3   -1   -11
// -2    0   -3    22

    int M = 3 ;
    int N = 4 ;
    double data[] = { 1, 2, -2,   2, 3, 0,   -1, -1, -3,   -4, -11, 22 } ;
    matprint( data, M, N ) ;
    rref( data, M, N ) ;
    matprint( data, M, N ) ;

    return 0 ;
}


void matprint( double *data, int M, int N ) {
    double *p = data ;
    for( int r=0 ; r<M ; r++ ) { 
        for( int c=0 ; c<N ; c++ ) {
            std::cout << p[c*M] << " " ;
        }
        std::cout << std::endl ;
        p++ ;
    }
    std::cout << std::endl ;
}