
#include "rref.hpp"

extern "C" {

    void rref( double *data, const size_t M, const size_t N ) {

        double *p = data ;

        int lead = 0 ;
    
        for( int rix=0; rix<M; rix++ ) {
            if( lead >= N )
                return ;
            int ix = rix;
            p = data + lead * M ;            
            while( *p == 0 ) {
                ix++; 
                p++ ;
                if( ix == M ) {
                    ix = rix;
                    lead++;
                    if( lead == N )
                        return;
                }
            }

            // MtxSwapRows(m, iix, rix );

            // MtxNormalizeRow(m, rix, lead );
            p = data + rix + lead*M ;
            double lead_value = *p ;
            for( int i=lead ; i<N ; i++ ) {
                *p /= lead_value ;
                p += M ;
            }

            for( int iix=0; iix<M; iix++) {
                if ( iix != rix ) {
                    p = data + iix + (lead * M) ;
                    // double lv = MtxGet(m, iix, lead );
                    double lead_value = *p ;
                    // MtxMulAndAddRows(m,iix, rix, -lv) ;
                    double *drow = data + iix ;
                    double *srow = data + rix ;
                    for( int i=0 ; i<N ; i++ ) {
                        *drow += lead_value * *srow ;
                        drow += M ;
                        srow += M ;
                    }
                }
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

}