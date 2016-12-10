/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Gemm_panmat( int nb_alg, int transa, int transb, 
		      PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
                      PLA_Obj beta,  PLA_Obj C )
{
  PLA_Obj 
    C_L = NULL, C_1 = NULL, C_1_dup = NULL, A_L = NULL,  A_1 = NULL, 
    A_1_dup = NULL, one = NULL, zero = NULL;
  int 
    size;

  PLA_Local_scal( beta, C );   

  PLA_Obj_global_width(  C, &size );
  nb_alg = min( size, nb_alg );

  PLA_Obj_view_all( C, &C_L );        PLA_Obj_view_all( A, &A_L );

  if ( transb == PLA_NO_TRANS ){
    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &A_1_dup );
    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &C_1_dup );

    PLA_Create_constants_conf_to( C, NULL, &zero, &one );

    while ( TRUE ){
      PLA_Obj_global_length( C_L, &size );
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      PLA_Obj_horz_split_2( C_L, size, &C_1, 
                                        &C_L );
      PLA_Obj_set_orientation( C_1, PLA_PROJ_ONTO_ROW );

      if ( transa == PLA_NO_TRANS ) {
	PLA_Obj_horz_split_2( A_L, size, &A_1, 
                                          &A_L );
	PLA_Obj_set_orientation( A_1, PLA_PROJ_ONTO_ROW );
      }
      else 
	PLA_Obj_vert_split_2( A_L, size, &A_1, &A_L );

      if ( size != nb_alg ) {
	PLA_Obj_vert_split_2( A_1_dup, size, &A_1_dup, PLA_DUMMY );
	PLA_Obj_horz_split_2( C_1_dup, size, &C_1_dup, 
                                              PLA_DUMMY );
      }

      PLA_Copy( A_1, A_1_dup );

      if ( transa == PLA_NO_TRANS ) 
	PLA_Local_gemm( PLA_TRANS, PLA_NO_TRANS, alpha, A_1_dup, B,
                                                    zero,  C_1_dup );
      else
	PLA_Local_gemm( transa, PLA_NO_TRANS,     alpha, A_1_dup, B,
                                                    zero,  C_1_dup );

      PLA_Reduce_x( PLA_SHAPE_GENERAL, C_1_dup, one, C_1 );   
    }
  }
  else 
    PLA_Abort( "transb == PLA_TRANS not yet implemented", __LINE__, __FILE__ );

  PLA_Obj_free( &C_L);   PLA_Obj_free( &C_1 );   PLA_Obj_free( &C_1_dup );
  PLA_Obj_free( &A_L);   PLA_Obj_free( &A_1 );   PLA_Obj_free( &A_1_dup );
  PLA_Obj_free( &one );  PLA_Obj_free( &zero );

  return PLA_SUCCESS;
}


