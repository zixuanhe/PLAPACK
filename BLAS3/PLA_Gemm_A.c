/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Gemm_A( int nb_alg, int transa, int transb, 
		      PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
                      PLA_Obj beta,  PLA_Obj C )
/*
  Purpose : Parallel matrix multiplication

  IN     nb_alg      integer, algorithmic block size
  IN     transa      integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     transb      integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     alpha       multiscalar, scale factor for A and B
  IN     A           matrix
  IN     B           matrix
  IN     beta        multiscalar, scale factor for C
  IN/OUT C           matrix, overwritten with 
                     C <- alpha * A  * B  + beta * C  or 
                     C <- alpha * A  * B' + beta * C  or 
                     C <- alpha * A' * B  + beta * C  or 
                     C <- alpha * A  * B' + beta * C

  This version assumes that A constitutes the bulk of the data.  Thus,
  the implementation leaves A in place, communicating B and C.
  
  For the case where C <- alpha * A * B + beta * C is to be computed,
  we partition C = ( C_0 ... C_(K-1) ) and B = ( B_0 ...  B_(k-1) )
  and iterate computing C_j = alpha * A * B_j + C_j.  In other words,
  the computation is set up as a sequence of matrix-panel( of columns) 
  multiplies.

  More precisely, the algorithm is given by

  ******************************************************************
      
      C <- beta * C
      Partition  B = ( B_F || B_L ) and C = ( C_F || C_L ) 
     	     where B_F is k x 0 and C_F is m x 0
      while C_L is not m x 0 
         determine block size b
         Partition ( B_F || B_L ) = ( B_0 || B_1 | B_2 ) 
                   where B_0 = B_F and B_1 has width b
              and  ( C_F || C_L ) = ( C_0 || C_1 | C_2 ) 
                   where C_0 = C_F and C_1 has width b
          Update C_1 <- alpha * A * B_1 + C_1    (matrix-panel mult.)
          Continue with
                   ( B_F || B_L ) = ( B_0 | B_1 || B_2 ) 
	      and  
                   ( C_F || C_L ) = ( C_0 | C_1 || C_2 ) 
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
    B_L = NULL,  B_1 = NULL, B_1_dup = NULL, 
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
  PLA_Obj_view_all( B, &B_L );

  if ( transa == PLA_NO_TRANS ){
    /* Create duplicated projected multivectors to hold copies of
       panels of B and local contributions to C_j */
    PLA_Pmvector_create_conf_to( 
          A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &C_1_dup );
    PLA_Pmvector_create_conf_to( 
          A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &B_1_dup );

    while ( TRUE ){
      PLA_Obj_global_width( C_L, &size );

      /* Last panel may not be of full width.  If zero, we are done */
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      /* Partition off next panel from C_L */
      PLA_Obj_vert_split_2( C_L, size, &C_1, &C_L );

      /* Partition off next panel from B_L */
      if ( transb == PLA_NO_TRANS ) 
	PLA_Obj_vert_split_2( B_L, size, &B_1, &B_L );
      else {
	PLA_Obj_horz_split_2( B_L, size, &B_1, 
                                          &B_L );
	PLA_Obj_set_orientation( B_1, PLA_PROJ_ONTO_ROW );
      }

      /* If we are on the last panel, it may be necessary to 
	 resize the duplicated projected multivectors. */
      if ( size != nb_alg ) {
	PLA_Obj_vert_split_2( C_1_dup, size, &C_1_dup, PLA_DUMMY );
	PLA_Obj_horz_split_2( B_1_dup, size, &B_1_dup, 
                                              PLA_DUMMY );
      }

      /* Duplicate the current panel B */
      PLA_Copy( B_1, B_1_dup );

      /* Compute local contributions to current panel of C */
      if ( transb == PLA_NO_TRANS ) 
	PLA_Local_gemm( PLA_NO_TRANS, PLA_TRANS, alpha, A, B_1_dup, 
                                                    zero,  C_1_dup );
      else
	PLA_Local_gemm( PLA_NO_TRANS, transb,     alpha, A, B_1_dup, 
                                                    zero,  C_1_dup );

      /* Add all contributions into current panel of C */
      PLA_Reduce_x( PLA_SHAPE_GENERAL, C_1_dup, one, C_1 );   
    }
  }
  else { /* transa != PLA_NO_TRANS */
    /* Create duplicated projected multivectors to hold copies of
       panels of B and local contributions to C_j */
    PLA_Pmvector_create_conf_to( 
          A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &B_1_dup );
    PLA_Pmvector_create_conf_to( 
          A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &C_1_dup );

    while ( TRUE ){
      PLA_Obj_global_width( C_L, &size );

      /* Last panel may not be of full width.  If zero, we are done */
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      /* Partition off next panel from C_L */
      PLA_Obj_vert_split_2( C_L, size, &C_1, &C_L );

      /* Partition off next panel from B_L */
      if ( transb == PLA_NO_TRANS ) 
	PLA_Obj_vert_split_2( B_L, size, &B_1, &B_L );
      else {
	PLA_Obj_horz_split_2( B_L, size, &B_1, 
                                          &B_L );
	PLA_Obj_set_orientation( B_1, PLA_PROJ_ONTO_ROW );
      }

      /* If we are on the last panel, it may be necessary to 
	 resize the duplicated projected multivectors. */
      if ( size != nb_alg ) {
	PLA_Obj_vert_split_2( B_1_dup, size, &B_1_dup, PLA_DUMMY );
	PLA_Obj_horz_split_2( C_1_dup, size, &C_1_dup, 
                                              PLA_DUMMY );
      }

      /* Duplicate the current panel B */
      PLA_Copy( B_1, B_1_dup );

      /* Compute local contributions to current panel of C */
      if ( transb == PLA_NO_TRANS ) 
	PLA_Local_gemm( PLA_TRANS, 
			 ( transa == PLA_TRANS ? PLA_NO_TRANS : PLA_CONJ ),
			 alpha, B_1_dup, A, zero, C_1_dup );
      else
	PLA_Local_gemm( transb, 
			 ( transa == PLA_TRANS ? PLA_NO_TRANS : PLA_CONJ ),
			 alpha, B_1_dup, A, zero,  C_1_dup );

      /* Add all contributions into current panel of C */
      PLA_Reduce_x( PLA_SHAPE_GENERAL, C_1_dup, one, C_1 );   
    }
  }
    
  /* free temporary objects and views */
  PLA_Obj_free( &C_L);   PLA_Obj_free( &C_1 );   PLA_Obj_free( &C_1_dup );
  PLA_Obj_free( &B_L);   PLA_Obj_free( &B_1 );   PLA_Obj_free( &B_1_dup );
  PLA_Obj_free( &one );  PLA_Obj_free( &zero );

  return PLA_SUCCESS;
}


