/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Gemm_panpan( int nb_alg, int transa, int transb, 
		      PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
                      PLA_Obj beta,  PLA_Obj C )
{
  PLA_Obj 
    A_L = NULL, A_1 = NULL, A_1_dup = NULL, B_L = NULL,  B_1 = NULL, 
    B_1_dup = NULL, one = NULL;
  int 
    size;

  PLA_Local_scal( beta, C );   
  PLA_Create_constants_conf_to( C, NULL, NULL, &one );

  if ( transa == PLA_NO_TRANS ) 
    PLA_Obj_global_width(  A, &size );
  else                          
    PLA_Obj_global_length( A, &size );
  nb_alg = min( size, nb_alg );

  PLA_Obj_view_all( A, &A_L );        PLA_Obj_view_all( B, &B_L );

  PLA_Pmvector_create_conf_to( 
          C, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &A_1_dup );
  PLA_Pmvector_create_conf_to( 
          C, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &B_1_dup );


  while ( TRUE ){
    if ( transa == PLA_NO_TRANS ) 
      PLA_Obj_global_width(  A_L, &size );
   else                          
      PLA_Obj_global_length( A_L, &size );
    if ( ( size = min( size, nb_alg ) ) == 0 ) break;

    if ( transa == PLA_NO_TRANS ) 
      PLA_Obj_vert_split_2( A_L, size, &A_1, &A_L );
    else {
      PLA_Obj_horz_split_2( A_L, size, &A_1, 
                                        &A_L );
      PLA_Obj_set_orientation( A_1, PLA_PROJ_ONTO_ROW );
    }

    if ( transb == PLA_NO_TRANS ) {
      PLA_Obj_horz_split_2( B_L, size, &B_1, 
                                        &B_L );
      PLA_Obj_set_orientation( B_1, PLA_PROJ_ONTO_ROW );
    }
    else 
      PLA_Obj_vert_split_2( B_L, size, &B_1, &B_L );

    if ( size != nb_alg ) {
      PLA_Obj_vert_split_2( A_1_dup, size, &A_1_dup, PLA_DUMMY );
      PLA_Obj_horz_split_2( B_1_dup, size, &B_1_dup, 
                                            PLA_DUMMY );
    }

    PLA_Copy( A_1, A_1_dup );       
    PLA_Copy( B_1, B_1_dup );

    if ( transa == PLA_CONJ_TRANS )   PLA_Conjugate( A_1_dup );
    if ( transb == PLA_CONJ_TRANS )   PLA_Conjugate( B_1_dup );

    PLA_Local_gemm( PLA_NO_TRANS, PLA_NO_TRANS, alpha, A_1_dup, B_1_dup, 
                                                   one,   C );
  }

  PLA_Obj_free( &A_L);   PLA_Obj_free( &A_1 );   PLA_Obj_free( &A_1_dup );
  PLA_Obj_free( &B_L);   PLA_Obj_free( &B_1 );   PLA_Obj_free( &B_1_dup );
  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}


