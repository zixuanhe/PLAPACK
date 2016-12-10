/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Gemm_C( int nb_alg, int transa, int transb, 
		      PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
                      PLA_Obj beta,  PLA_Obj C )
/*
int PLA_Gemm_C ( int nb_alg
                 int transa,
                 int transb,
                 PLA_Obj alpha,
                 PLA_Obj A,
                 PLA_Obj B,
                 PLA_Obj beta,
                 PLA_Obj C)

  Purpose : Parallel matrix multiplication

  IN     nb_alg      integer, algorithmic block size to be used
  IN     transa      integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     transb      integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     alpha       multiscalar, scale factor for A and B
  IN     A           matrix, factor in product
  IN     B           matrix, factor in product
  IN     beta        multiscalar, scale factor for C
  IN/OUT C           matrix, overwritten with 
                     C <- alpha * A  * B  + beta * C  or 
                     C <- alpha * A  * B' + beta * C  or 
                     C <- alpha * A' * B  + beta * C  or 
                     C <- alpha * A  * B' + beta * C
  		     
  This version assumes that C constitutes the bulk of the data.  Thus,
  the implementation leaves C in place, communicating A and B.
  
  For the case where C <- alpha * A * B + beta * C is to be computed,
  we partition A = ( A_0 ... A_(K-1) ) and B = /   B_0   \
                                               |    :    |
                                               \ B_(k-1) /
  and iterate updating C = alpha * A_j * B_j + C.  In other words,
  the computation is set up as a sequence of rank-k updates, with 
  k = nb_alg.

  More precisely, the algorithm is given by

  ******************************************************************
      
      C <- beta * C
      Partition  A = ( A_F || A_L ) and B = / B_F \
                                            | === |
                                            \ B_L /
     	     where A_F is m x 0 and B_F is 0 x n
      while A_L is not m x 0 
         determine block size b
         Partition ( A_F || A_L ) = ( A_0 || A_1 | A_2 ) 
                   where A_0 = A_F and A_1 has width b
              and  / B_F \    / B_0 \
                   | === | =  | === |
                   \ B_L /    | B_1 |
                              | --- | 
                              \ B_2 /
                   where B_0 = B_F and B_1 has length b
          Update C <- alpha * A_1 * B_1 + C    (rank-b update)
          Continue with
                   ( A_F || A_L ) = ( A_0 | A_1 || A_2 ) 
    
              and  / B_F \    / B_0 \
                   | === | =  | --- |
                   \ B_L /    | B_1 |
                              | === | 
                              \ B_2 /
      endwhile
              
  ******************************************************************

  Appropriate changes need to be made depending on transa and transb

  NOTE:  For details on how to implement matrix-matrix multiplication,
         see
	 
	 R. van de Geijn, Using PLAPACK, The MIT Press, 1997.

	 R. van de Geijn and J. Watts, 
	 "SUMMA: Scalable Universal Matrix Multiplication Algorith,"
	 Concurrency: Practice and Experience, Vol 9 (4), pp. 255-274
	 (April 1997).
          	 
	 G. Morrow, J. Gunnels, C. Lin, and R. van de Geijn,
	 "A Flexible Class of Parallel Matrix Multiplication Algorithms,"
	 Proceedings of IPPS98, pp. 110-116, 1998.
*/
{
  PLA_Obj 
    A_L = NULL, A_1 = NULL, A_1_dup = NULL, 
    B_L = NULL, B_1 = NULL, 
    B_1_dup = NULL, one = NULL;
  int 
    size;

  /* Scale C <- beta * C */
  PLA_Local_scal( beta, C );   
  
  /* Create usual constants */
  PLA_Create_constants_conf_to( C, NULL, NULL, &one );

  /* Create duplicated projected multivectors to hold copies of
     panels of A and B */
  if ( transa == PLA_NO_TRANS ) 
    PLA_Obj_global_width(  A, &size );
  else                          
    PLA_Obj_global_length( A, &size );
  nb_alg = min( size, nb_alg );

  PLA_Pmvector_create_conf_to( 
          C, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &A_1_dup );
  PLA_Pmvector_create_conf_to( 
          C, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &B_1_dup );

  /* Let A_L and B_L view the parts of A and B yet to be used
     (here "L" stands for Last) */
  PLA_Obj_view_all( A, &A_L );        
  PLA_Obj_view_all( B, &B_L );

  /* Loop over A_L and B_L, peeling off a next panel from each, and
     updating C with those panels. */

  while ( TRUE ){
    if ( transa == PLA_NO_TRANS ) 
      /* Check how many columns are left in A_L */
      PLA_Obj_global_width(  A_L, &size );
   else                          
      /* Check how many rows are left in A_L */
      PLA_Obj_global_length( A_L, &size );

    /* Last panel may not be of full width.  If zero, we are done */
    if ( ( size = min( size, nb_alg ) ) == 0 ) break;

    /* Partition off next panel from A_L */
    if ( transa == PLA_NO_TRANS ) 
      /* Let A_L = ( A_1 | A_L ) */
      PLA_Obj_vert_split_2( A_L, size, &A_1, &A_L );
    else {
      /* Let A_L = / A_1 \
                   \ A_L / */
      PLA_Obj_horz_split_2( A_L, size, &A_1, 
                                        &A_L );
      PLA_Obj_set_orientation( A_1, PLA_PROJ_ONTO_ROW );
    }

    /* Partition off next panel from B_L */
    if ( transb == PLA_NO_TRANS ) {
      /* Let B_L = / B_1 \
                   \ B_L / */
      PLA_Obj_horz_split_2( B_L, size, &B_1, 
                                        &B_L );
      PLA_Obj_set_orientation( B_1, PLA_PROJ_ONTO_ROW );
    }
    else 
      /* Let B_L = ( B_1 | B_L ) */
      PLA_Obj_vert_split_2( B_L, size, &B_1, &B_L );

    /* If we are on the last panel, it may be necessary to 
       resize the duplicated projected multivectors. */
    if ( size != nb_alg ) {
      PLA_Obj_vert_split_2( A_1_dup, size, &A_1_dup, PLA_DUMMY );
      PLA_Obj_horz_split_2( B_1_dup, size, &B_1_dup, 
                                            PLA_DUMMY );
    }

    /* Duplicate the current panels of A and B */
    PLA_Copy( A_1, A_1_dup );       
    PLA_Copy( B_1, B_1_dup );

    /* Since copy may transpose, it may be necessary to conjugate */
    if ( transa == PLA_CONJ_TRANS )   PLA_Conjugate( A_1_dup );
    if ( transb == PLA_CONJ_TRANS )   PLA_Conjugate( B_1_dup );

    /* Perform local part of update C <- alpha * A_1 * B_1 + C */
    PLA_Local_gemm( PLA_NO_TRANS, PLA_NO_TRANS, alpha, A_1_dup, B_1_dup, 
                                                   one,   C );
  }

  /* free temporary objects and views */
  PLA_Obj_free( &A_L);   PLA_Obj_free( &A_1 );   PLA_Obj_free( &A_1_dup );
  PLA_Obj_free( &B_L);   PLA_Obj_free( &B_1 );   PLA_Obj_free( &B_1_dup );
  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}


