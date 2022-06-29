#include <iostream>
#include <string.h>
#include <iomanip>

#include "rref.hpp"

extern "C" {
//https://rosettacode.org/wiki/Reduced_row_echelon_form#C

    void rref( double *data, const size_t M, const size_t N ) {

        int lead = 0 ;
    
        for( int rix=0; rix<M; rix++ ) {
            if( lead >= N ) return ;
            int iix = rix ;

            while( abs( data[iix + lead*M ] ) < 1e-11 ) {
                iix++ ; 
                if( iix == M ) {
                    iix = rix ;
                    lead++ ;
                    if( lead == N ) return;
                }
            }

            // MtxSwapRows(m, iix, rix );
            if( iix != rix ) {
                double *a = data + rix ;
                double *b = data + iix ;
                for( int i=0 ; i<N ; i++ ) {
                    double tmp = *a ;
                    *a = *b ;
                    *b = tmp ;
                    a += M ;
                    b += M ;
                }
            }

            // MtxNormalizeRow(m, rix, lead );
            double *p = data + rix + lead*M ;
            double lead_value = *p ;
            if( abs(lead_value) > 1e-12 ) {
                p = data + rix ;
                for( int i=0 ; i<N ; i++ ) {
                    *p /= lead_value ;
                    p += M ;
                }
            }
            // matprint( data, M, N ) ;

            for( int iix=0; iix<M; iix++) {
                if ( iix != rix ) {
                    p = data + iix + (lead * M) ;
                    // double lv = MtxGet(m, iix, lead );
                    double lead_value = -(*p) ;
                    // MtxMulAndAddRows(m,iix, rix, -lv) ;
                    double *drow = data + iix ;
                    double *srow = data + rix ;
                    for( int i=0 ; i<N ; i++ ) {
                        *drow += lead_value * *srow ;
                        drow += M ;
                        srow += M ;
                    }
                }
                // matprint( data, M, N ) ;
            }
            lead++;
        }
    }

// void MtxMulAndAddRows(Matrix m, int ixrdest, int ixrsrc, EL_Type mplr)
// {
//     int ix;
//     EL_Type *drow, *srow;
//     drow = m->mtx[ixrdest];
//     srow = m->mtx[ixrsrc];
//     for (ix=0; ix<m->dim_x; ix++) 
//         drow[ix] += mplr * srow[ix];
// //	printf("Mul row %d by %d and add to row %d\n", ixrsrc, mplr, ixrdest);
// //	MtxDisplay(m);
// }

/**
 * @ref https://adamdhalla.medium.com/linear-algebra-3-solving-ax-0-free-variables-and-pivot-variables-and-echelon-form-ce8cca411e60
 * 
 * @brief 
 * 
 * @param data 
 * @param M 
 * @param N 
 */

    int nullspace( double *data, const size_t M, const size_t N ) {
        rref( data, M, N ) ;
        bool pivots[N] ;

        int numPivots = 0 ;        
        int lead = 0 ;
        double *p = data ;
        for( int i=0 ; i<N && i<M; i++ ) {
            if( abs(*p) > 1e-11 ) {
                pivots[i] = true ;
                numPivots++ ;
                p += M+1 ; // skip down diagonal
                lead++ ;
            } else {
                pivots[i] = false ;
                while( abs(*p) < 1e-11 ) {
                    p += M ;
                    lead++ ;
                    if( lead >= N ) goto escape ;
                }
            }
        }
        escape:
        int numFreeVariables = N - numPivots ;
        double nullspace[N*numFreeVariables] ;
        memset( nullspace, 0, sizeof( nullspace) ) ;
        p = nullspace ;
        for( int i=0 ; i<N ; i++ ) {
            if( !pivots[i] ) { 
                *p = 1 ;
                p += N ;
            }
            p++ ;
        }

        // matprint( nullspace, N, numFreeVariables ) ;

        int pivot = 0 ;
        p = data ;
        for( int i=0 ; i<N ; i++ ) {  // down nullspace
            if( pivots[i] ) { 
                for( int f=0 ; f<numFreeVariables ; f++ ) { // across nullspace
                    double total = 0 ;
                    double *p2 = nullspace + f*N ;
                    for( int j=0 ; j<N ; j++ ) {  // across rref (data)
                        total += p[M*j] * -(*p2) ;
                        p2++ ;
                    }                
                    nullspace[i+f*N] = total ;
                }
                pivot++ ;
                p++ ;
                if( pivot >= numPivots ) break ;
            }
            // matprint( nullspace, N, numFreeVariables ) ;
        }
        // matprint( nullspace, N, numFreeVariables ) ;
        memcpy( data, nullspace, sizeof(nullspace) ) ;
        return numFreeVariables ;
    }


    int colspace( double *data, const size_t M, const size_t N ) {
        double tmp[M*N] ;
        memcpy( tmp, data, sizeof(tmp) ) ;

        rref( data, M, N ) ;
        int pivots[N] ;

        int numPivots = 0 ;        
        int lead = 0 ;
        double *p = data ;
        for( int i=0 ; i<N && i<M; i++ ) {
            if( abs(*p) > 1e-11 ) {
                pivots[numPivots] = lead ;
                numPivots++ ;
                p += M+1 ; // skip down diagonal
                lead++ ;
            } else {
                while( abs(*p) < 1e-11 ) {
                    p += M ;
                    lead++ ;
                    if( lead >= N ) goto escape ;
                }
            }
        }
        escape:
        
        double *d = data ;
        for( int i=0 ; i<numPivots ; i++ ) {
            double *s = tmp + M*pivots[i] ;
            for( int j=0 ; j<M ; j++ ) {
                *d++ = *s++ ;
            }
        }
        return numPivots ;
    }


    void matprint( double *data, int M, int N ) {
        double *p = data ;
        for( int r=0 ; r<M ; r++ ) { 
            for( int c=0 ; c<N ; c++ ) {
                std::cout << std::setw(5) << p[c*M] << " " ;
            }
            std::cout << std::endl ;
            p++ ;
        }
        std::cout << std::endl ;
    }
}