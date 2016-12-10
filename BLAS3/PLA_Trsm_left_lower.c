/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Trsm_left_lower( int transa, int diag,
			  PLA_Obj alpha, PLA_Obj A, PLA_Obj B )
{
  PLA_Obj 
    A_BR = NULL, A_11 = NULL, A_21 = NULL, A_21_dup = NULL, A_11_dup = NULL, 
    A_TL = NULL,              A_10 = NULL, A_10_dup = NULL, 
    B_B = NULL,  B_T = NULL, B_1 = NULL, C_1_dup = NULL,
    minus_one = NULL, one = NULL;

  int 
    nb_alg, size;

  PLA_Template
    templ;

  PLA_Local_scal( alpha, B );   
  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );
  
  PLA_Obj_template( A, &templ );

  if ( transa == PLA_NO_TRANS ){
    PLA_Obj_global_length( B, &size );
    PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );
    nb_alg = min( size, nb_alg );

    PLA_Obj_view_all( A, &A_BR );        PLA_Obj_view_all( B, &B_B );

    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &C_1_dup );

    while ( TRUE ){
      PLA_Obj_global_length( A_BR, &size );
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      PLA_Obj_split_4( A_BR, size, size, &A_11, PLA_DUMMY,
                                          &A_21, &A_BR );

      PLA_Obj_horz_split_2( B_B, size,   &B_1,
			                  &B_B );

      PLA_Obj_set_orientation( B_1, PLA_PROJ_ONTO_ROW );

      PLA_Mscalar_create_conf_to( A_11, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &A_11_dup );

      PLA_Pmvector_create_conf_to( B_B, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
				    size, &A_21_dup );

      if ( size != nb_alg ){
	PLA_Obj_horz_split_2( C_1_dup,   size, &C_1_dup, 
                                                PLA_DUMMY );
      }

      PLA_Copy( B_1, C_1_dup );       
      PLA_Copy( A_11, A_11_dup );
      PLA_Copy( A_21, A_21_dup );

      PLA_Local_trsm( PLA_SIDE_LEFT, PLA_LOWER_TRIANGULAR,
		       PLA_NO_TRANSPOSE, diag,
		       one, A_11_dup, C_1_dup );
      
      PLA_Copy( C_1_dup, B_1 );
      
      PLA_Local_gemm( PLA_NO_TRANS, PLA_NO_TRANS,
		       minus_one, A_21_dup, C_1_dup, one, B_B );
    }
  }
  else{ /* transa != PLA_NO_TRANSPOSE */
    PLA_Obj_global_length( B, &size );
    PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );
    nb_alg = min( size, nb_alg );

    PLA_Obj_view_all( A, &A_TL );        PLA_Obj_view_all( B, &B_T );

    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &C_1_dup );

    while ( TRUE ){
      PLA_Obj_global_length( A_TL, &size );
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      PLA_Obj_split_4( A_TL, -size, -size, &A_TL, PLA_DUMMY,
                                            &A_10, &A_11 );

      PLA_Obj_set_orientation( A_10, PLA_PROJ_ONTO_ROW );

      PLA_Obj_horz_split_2( B_T, -size,   &B_T,
			                   &B_1 );

      PLA_Obj_set_orientation( B_1, PLA_PROJ_ONTO_ROW );

      PLA_Mscalar_create_conf_to( A_11, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &A_11_dup );

      PLA_Pmvector_create_conf_to( B_T, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
				    size, &A_10_dup );

      if ( size != nb_alg ){
	PLA_Obj_horz_split_2( C_1_dup,   size, &C_1_dup, 
                                                PLA_DUMMY );
      }

      PLA_Copy( B_1, C_1_dup );       
      PLA_Copy( A_11, A_11_dup );
      PLA_Copy( A_10, A_10_dup );

      PLA_Local_trsm( PLA_SIDE_LEFT, PLA_LOWER_TRIANGULAR,
		       transa, diag,
		       one, A_11_dup, C_1_dup );
      
      PLA_Copy( C_1_dup, B_1 );
      
      if ( transa == PLA_CONJ_TRANS )
	PLA_Conjugate( A_10_dup );
      PLA_Local_gemm( PLA_NO_TRANS, PLA_NO_TRANS,
		       minus_one, A_10_dup, C_1_dup, one, B_T );
    }
  }

  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_TL );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_21 );
  PLA_Obj_free( &A_10 );
  PLA_Obj_free( &A_11_dup );
  PLA_Obj_free( &A_21_dup );
  PLA_Obj_free( &A_10_dup );
  PLA_Obj_free( &B_B );
  PLA_Obj_free( &B_T );
  PLA_Obj_free( &B_1 );
  PLA_Obj_free( &C_1_dup );
  PLA_Obj_free( &one );
  PLA_Obj_free( &minus_one );


  return PLA_SUCCESS;
}


