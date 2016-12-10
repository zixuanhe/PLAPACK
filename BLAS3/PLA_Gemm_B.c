/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Gemm_B( int nb_alg, int transa, int transb, 
		 PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
                 PLA_Obj beta,  PLA_Obj C )
/*
int PLA_Gemm_B ( int nb_alg,
                 int transa,
                 int transb,
                 PLA_Obj alpha,
                 PLA_Obj A,
                 PLA_Obj B,
                 PLA_Obj beta,
                 PLA_Obj C)

  Purpose : Parallel matrix multiplication

  IN     nb_alg      integer, algorithmic block size
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

  This version assumes that B constitutes the bulk of the data.  Thus,
  the implementation leaves B in place, communicating A and C.
  
  For the case where C <- alpha * A * B + beta * C is to be computed,
  we partition C = /   C_0   \ and /   A_0   \
                   |    :    |     |    :    |
                   \ C_(K-1) /     \ A_(k-1) /
  and iterate computing C_i = alpha * A_i * B + C_i.  In other words,
  the computation is set up as a sequence of panel( of rows)-matrix
  multiplies.

  More precisely, the algorithm is given by

  ******************************************************************
      
      C <- beta * C
      Partition  A = / A_F \ and C = / C_F \
                     | === |         | === |
                     \ A_L /         \ C_L /
     	     where A_F is 0 x k and C_F is 0 x n
      while A_L is not 0 x k 
         determine block size b
         Partition 
              / A_F \    / A_0 \       / C_F \    / C_0 \
              | === | =  | === |  and  | === | =  | === |
              \ A_L /    | A_1 |       \ C_L /    | C_1 |
                         | --- |                  | --- | 
                         \ A_2 /                  \ C_2 /
                   where A_0 = A_F and A_1 has length b
                         C_0 = C_F and C_1 has length b
          Update C_1 <- alpha * A_1 * B + C_1    (panel-matrix mult.)
          Continue with
              / A_F \    / A_0 \       / C_F \    / C_0 \
              | === | =  | --- |  and  | === | =  | --- |
              \ A_L /    | A_1 |       \ C_L /    | C_1 |
                         | === |                  | === | 
                         \ A_2 /                  \ C_2 /
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
    C_L = NULL, C_1 = NULL, C_1_dup = NULL, 
    A_L = NULL,  A_1 = NULL, A_1_dup = NULL, 
    one = NULL, zero = NULL;
  int 
    size;

  /* Scale C <- beta * C */
  PLA_Local_scal( beta, C );   

  /* Create usual constants */
  PLA_Create_constants_conf_to( C, NULL, &zero, &one );

  PLA_Obj_global_width(  C, &size );
  nb_alg = min( size, nb_alg );

  /* Let B_L and C_L view the parts of B and C yet to be used/updated
     (here "L" stands for Last) */
  PLA_Obj_view_all( C, &C_L ); 
  PLA_Obj_view_all( A, &A_L );

  if ( transb == PLA_NO_TRANS ){
    /* Create duplicated projected multivectors to hold copies of
       panels of A and local contributions to C_i */
    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &A_1_dup );
    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &C_1_dup );

    while ( TRUE ){
      PLA_Obj_global_length( C_L, &size );

      /* Last panel may not be of full width.  If zero, we are done */
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      /* Partition off next panel from C_L */
      PLA_Obj_horz_split_2( C_L, size, &C_1, 
                                       &C_L );
      PLA_Obj_set_orientation( C_1, PLA_PROJ_ONTO_ROW );

      /* Partition off next panel from A_L */
      if ( transa == PLA_NO_TRANS ) {
	PLA_Obj_horz_split_2( A_L, size, &A_1, 
                                         &A_L );
	PLA_Obj_set_orientation( A_1, PLA_PROJ_ONTO_ROW );
      }
      else 
	PLA_Obj_vert_split_2( A_L, size, &A_1, &A_L );

      /* If we are on the last panel, it may be necessary to 
	 resize the duplicated projected multivectors. */
      if ( size != nb_alg ) {
	PLA_Obj_vert_split_2( A_1_dup, size, &A_1_dup, PLA_DUMMY );
	PLA_Obj_horz_split_2( C_1_dup, size, &C_1_dup, 
                                              PLA_DUMMY );
	PLA_Obj_set_orientation( C_1_dup, PLA_PROJ_ONTO_ROW );
      }

      /* Duplicate the current panel A */
      PLA_Copy( A_1, A_1_dup );

      /* Compute local contributions to current panel of C */
      if ( transa == PLA_NO_TRANS ) 
	PLA_Local_gemm( PLA_TRANS, PLA_NO_TRANS, alpha, A_1_dup, B,
                                                    zero,  C_1_dup );
      else
	PLA_Local_gemm( transa, PLA_NO_TRANS,     alpha, A_1_dup, B,
                                                  zero,  C_1_dup );

      /* Add all contributions into current panel of C */
      PLA_Reduce_x( PLA_SHAPE_GENERAL, C_1_dup, one, C_1 );
    }
  }
  else {
    /* Create duplicated projected multivectors to hold copies of
       panels of A and local contributions to C_i */
    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &C_1_dup );
    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &A_1_dup );

    if ( transb == PLA_CONJ_TRANS ) PLA_Conjugate( B );

    while ( TRUE ){
      PLA_Obj_global_length( C_L, &size );

      /* Last panel may not be of full width.  If zero, we are done */
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      /* Partition off next panel from C_L */
      PLA_Obj_horz_split_2( C_L, size, &C_1, 
                                        &C_L );
      PLA_Obj_set_orientation( C_1, PLA_PROJ_ONTO_ROW );

      /* Partition off next panel from A_L */
      if ( transa == PLA_NO_TRANS ) {
	PLA_Obj_horz_split_2( A_L, size, &A_1, 
                                          &A_L );
	PLA_Obj_set_orientation( A_1, PLA_PROJ_ONTO_ROW );
      }
      else 
	PLA_Obj_vert_split_2( A_L, size, &A_1, &A_L );

      /* If we are on the last panel, it may be necessary to 
	 resize the duplicated projected multivectors. */
      if ( size != nb_alg ) {
	PLA_Obj_vert_split_2( C_1_dup, size, &C_1_dup, PLA_DUMMY );
	PLA_Obj_horz_split_2( A_1_dup, size, &A_1_dup, 
                                              PLA_DUMMY );
      }

      /* Duplicate the current panel A */
      PLA_Copy( A_1, A_1_dup );

      /* Compute local contributions to current panel of C */
      if ( transa == PLA_NO_TRANS )
	PLA_Local_gemm( PLA_NO_TRANS, PLA_TRANS, alpha, B, A_1_dup, 
                                                    zero,  C_1_dup );
      else 
	PLA_Local_gemm( PLA_NO_TRANS, transa, alpha, B, A_1_dup,
                                                       zero,  C_1_dup );

      /* Add all contributions into current panel of C */
      PLA_Reduce_x( PLA_SHAPE_GENERAL, C_1_dup, one, C_1 );   
    }

    if ( transb == PLA_CONJ_TRANS ) PLA_Conjugate( B );
  }

  /* free temporary objects and views */
  PLA_Obj_free( &C_L);   PLA_Obj_free( &C_1 );   PLA_Obj_free( &C_1_dup );
  PLA_Obj_free( &A_L);   PLA_Obj_free( &A_1 );   PLA_Obj_free( &A_1_dup );
  PLA_Obj_free( &one );  PLA_Obj_free( &zero );

  return PLA_SUCCESS;
}


