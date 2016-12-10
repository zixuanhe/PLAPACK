/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Gemm_matpan( int nb_alg, int transa, int transb, 
		      PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
                      PLA_Obj beta,  PLA_Obj C )
{
  PLA_Obj 
    C_L = NULL, C_1 = NULL, C_1_dup = NULL, B_L = NULL,  B_1 = NULL, 
    B_1_dup = NULL, one = NULL, zero = NULL;
  int 
    size;

  PLA_Local_scal( beta, C );   

  PLA_Obj_global_width(  C, &size );
  nb_alg = min( size, nb_alg );

  PLA_Obj_view_all( C, &C_L );        PLA_Obj_view_all( B, &B_L );

  PLA_Create_constants_conf_to( C, NULL, &zero, &one );

  if ( transa == PLA_NO_TRANS ){
    PLA_Pmvector_create_conf_to( 
          A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &C_1_dup );
    PLA_Pmvector_create_conf_to( 
          A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &B_1_dup );

    while ( TRUE ){
      PLA_Obj_global_width( C_L, &size );
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      PLA_Obj_vert_split_2( C_L, size, &C_1, &C_L );

      if ( transb == PLA_NO_TRANS ) 
	PLA_Obj_vert_split_2( B_L, size, &B_1, &B_L );
      else {
	PLA_Obj_horz_split_2( B_L, size, &B_1, 
                                          &B_L );
	PLA_Obj_set_orientation( B_1, PLA_PROJ_ONTO_ROW );
      }

      if ( size != nb_alg ) {
	PLA_Obj_vert_split_2( C_1_dup, size, &C_1_dup, PLA_DUMMY );
	PLA_Obj_horz_split_2( B_1_dup, size, &B_1_dup, 
                                              PLA_DUMMY );
      }

      PLA_Copy( B_1, B_1_dup );

      if ( transb == PLA_NO_TRANS ) 
	PLA_Local_gemm( PLA_NO_TRANS, PLA_TRANS, alpha, A, B_1_dup, 
                                                    zero,  C_1_dup );
      else
	PLA_Local_gemm( PLA_NO_TRANS, transb,     alpha, A, B_1_dup, 
                                                    zero,  C_1_dup );

      PLA_Reduce_x( PLA_SHAPE_GENERAL, C_1_dup, one, C_1 );   
    }
  }
  else { /* transa != PLA_NO_TRANS */
    PLA_Pmvector_create_conf_to( 
          A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &B_1_dup );
    PLA_Pmvector_create_conf_to( 
          A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &C_1_dup );

    while ( TRUE ){
      PLA_Obj_global_width( C_L, &size );
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      PLA_Obj_vert_split_2( C_L, size, &C_1, &C_L );

      if ( transb == PLA_NO_TRANS ) 
	PLA_Obj_vert_split_2( B_L, size, &B_1, &B_L );
      else {
	PLA_Obj_horz_split_2( B_L, size, &B_1, 
                                          &B_L );
	PLA_Obj_set_orientation( B_1, PLA_PROJ_ONTO_ROW );
      }

      if ( size != nb_alg ) {
	PLA_Obj_vert_split_2( B_1_dup, size, &B_1_dup, PLA_DUMMY );
	PLA_Obj_horz_split_2( C_1_dup, size, &C_1_dup, 
                                              PLA_DUMMY );
      }

      PLA_Copy( B_1, B_1_dup );

      if ( transb == PLA_NO_TRANS ) 
	PLA_Local_gemm( PLA_TRANS, 
			 ( transa == PLA_TRANS ? PLA_NO_TRANS : PLA_CONJ ),
			 alpha, B_1_dup, A, zero, C_1_dup );
      else
	PLA_Local_gemm( transb, 
			 ( transa == PLA_TRANS ? PLA_NO_TRANS : PLA_CONJ ),
			 alpha, B_1_dup, A, zero,  C_1_dup );

      PLA_Reduce_x( PLA_SHAPE_GENERAL, C_1_dup, one, C_1 );   
    }
  }
    
    

  PLA_Obj_free( &C_L);   PLA_Obj_free( &C_1 );   PLA_Obj_free( &C_1_dup );
  PLA_Obj_free( &B_L);   PLA_Obj_free( &B_1 );   PLA_Obj_free( &B_1_dup );
  PLA_Obj_free( &one );  PLA_Obj_free( &zero );

  return PLA_SUCCESS;
}


