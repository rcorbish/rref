
#pragma once 

#include <stdlib.h>

extern "C" {
    void rref( double *data, const size_t M, const size_t N ) ;
    int nullspace( double *data, const size_t M, const size_t N ) ;
    void matprint( double *data, int M, int N ) ;
}
