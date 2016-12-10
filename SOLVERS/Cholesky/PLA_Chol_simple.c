/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Chol_simple( int nb_alg, PLA_Obj A )
/*
  Purpose : Parallel Cholesky Factorization.  This particular version 
            assumes only the lower triangular portion of A is stored.
	    It only uses PLAPACK parallel BLAS calls, to illustrate how
            a simple implementation reflects the algorithm perfectly.

  IN     nb_alg      integer, algorithmic block size to be used
  IN/OUT A           matrix to be factored

  Algorith used:

  ******************************************************************
      
      Partition  A = / A_TL ||   *  \
                     | =====  ===== |
                     \ A_BL || A_BR /
     	     where A_TL is 0 x 0 
      while A_BR is not 0 x 0 
         Determine block size b
	 Partition  A = / A_TL ||   *  \    / A_00 ||  *   |  *   \ 
                        | =====  ===== | =  | ====    ====   ==== |
                        \ A_BL || A_BR /    | A_10 || A_11 |  *   |
                                            | ----    ----   ---- |
                                            \ A_20 || A_21 | A_22 /
     	        where A_00 = A_TL and A_11 is b x b 
         Update A_11 <- L_11 = Chol. Fact.( A_11 )
         Update A_21 <- L_21 = A_21 inv( L_11' )
	 Update A_22 <- A_22 - L_21 * L_21'
         Continue with
	            A = / A_TL ||   *  \    / A_00 |  *   ||  *   \ 
                        | =====  ===== | =  | ----   ----    ---- |
                        \ A_BL || A_BR /    | A_10 | A_11 ||  *   |
                                            | ====   ====    ==== |
                                            \ A_20 | A_21 || A_22 /
      endwhile
              
  ******************************************************************

  NOTE:  For details on how to implement parallel Cholesky factorization
         see
	 
	 R. van de Geijn, Using PLAPACK, The MIT Press, 1997.

	 G. Morrow and R. van de Geijn, "Zen and the Art of High-Performance
                 Parallel Computing," http://www.cs.utexas.edu/users/plapack
*/
{
  int       
    value = 0,
    size_top, size_left,
    owner_top, owner_left,
    size;

  PLA_Obj  
    ABR             = NULL,     
    A11             = NULL,     A21             = NULL,     
    one             = NULL,     minus_one       = NULL;

  /* Create usual duplicated scalar constants */
  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  /* View ABR = A */
  PLA_Obj_view_all( A, &ABR );

  while ( TRUE ) {
    /* Determine size of current panel.  Notice that we limit the size so that 
       A_11 resides on a single processor */
    PLA_Obj_split_size( ABR, PLA_SIDE_TOP,  &size_top,  &owner_top );
    PLA_Obj_split_size( ABR, PLA_SIDE_LEFT, &size_left, &owner_left );
    if ( 0 == ( size = min( min( size_top, size_left), nb_alg ) ) ) break;

    /* Partition A_BR = / A_11   *   \
                        \ A_21  A_22 / where A_11 is b x b.
       Notice that A_22 becomes A_BR in the next iteration, so we
       already view this as ABR. */
    PLA_Obj_split_4( ABR, size, size,   &A11, PLA_DUMMY,
                                         &A21, &ABR );

    /* Update A_11 <- L_11 = Chol. Fact.( A_11 ) */
    PLA_Local_chol( PLA_LOWER_TRIANGULAR, A11 );

    /* Update A_21 <- L_21 = A_21 inv( L_11' ) */
    PLA_Trsm( PLA_SIDE_RIGHT, PLA_LOWER_TRIANGULAR,
	       PLA_TRANSPOSE,  PLA_NONUNIT_DIAG,
	       one, A11, A21 );

    /* Update A_22 <- A_22 - L_21 * L_21' */
    PLA_Syrk( PLA_LOWER_TRIANGULAR, minus_one, A21, one, ABR );
  }

  /* free temporary objects and views */
  PLA_Obj_free( &ABR );        
  PLA_Obj_free( &A11 );        PLA_Obj_free( &A21 ); 
  PLA_Obj_free( &one );        PLA_Obj_free( &minus_one );
  
  return value;
}



