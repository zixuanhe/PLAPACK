/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Syrk_panpan( int nb_alg, int uplo, int transa,
		      PLA_Obj alpha, PLA_Obj A, 
                      PLA_Obj beta,  PLA_Obj C )
{
  PLA_Obj 
    A_L = NULL, A_1 = NULL, A_1_dup = NULL, A_1_trans_dup = NULL,
    one = NULL;
  int 
    size;

  PLA_Local_scal( beta, C );   
  PLA_Create_constants_conf_to( C, NULL, NULL, &one );

  if ( transa == PLA_NO_TRANS ) 
    PLA_Obj_global_width(  A, &size );
  else                          
    PLA_Obj_global_length( A, &size );
  nb_alg = min( size, nb_alg );

  PLA_Obj_view_all( A, &A_L );

  PLA_Pmvector_create_conf_to( 
          C, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &A_1_dup );
  PLA_Pmvector_create_conf_to( 
          C, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &A_1_trans_dup );


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

    if ( size != nb_alg ) {
      PLA_Obj_vert_split_2( A_1_dup,       size, &A_1_dup, PLA_DUMMY );
      PLA_Obj_horz_split_2( A_1_trans_dup, size, &A_1_trans_dup, 
                                                   PLA_DUMMY );
    }

    PLA_Copy( A_1,     A_1_dup );       
    PLA_Copy( A_1_dup, A_1_trans_dup );

    if ( transa == PLA_CONJ_TRANS )   PLA_Conjugate( A_1_dup );

    PLA_Syrk_perform_local_part( uplo, alpha, A_1_dup, A_1_trans_dup, one,   C );
  }

  PLA_Obj_free( &A_L);           PLA_Obj_free( &A_1 );   
  PLA_Obj_free( &A_1_dup );      PLA_Obj_free( &A_1_trans_dup );
  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}


