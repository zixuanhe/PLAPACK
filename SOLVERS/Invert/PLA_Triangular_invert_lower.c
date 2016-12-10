/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Triangular_invert_lower( PLA_Obj A )
{
  int 
    nb_alg, size, value = PLA_SUCCESS;

  PLA_Obj
    A_BR = NULL, A_11 = NULL, A_21 = NULL, A_21_dup = NULL, A_11_msc = NULL,
    Ainv = NULL, Ainv_BL = NULL, Ainv_10_11 = NULL, Ainv_10_11_dup = NULL,
    minus_one = NULL, zero = NULL, one = NULL;

  PLA_Template
    templ;

  PLA_Obj_template( A, &templ );
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  PLA_Create_constants_conf_to( A, &minus_one, &zero, &one );
  PLA_Obj_vert_split_2( A, 0, &Ainv_BL, &A_BR );

  
  while( TRUE ){
    PLA_Obj_global_length( A_BR, &size );
    if ( 0 == ( size = min( size, nb_alg ) ) ) break;

    PLA_Obj_split_4( A_BR, size, size, &A_11, PLA_DUMMY,
                                        &A_21, &A_BR );

    PLA_Obj_view_shift( Ainv_BL,       0,
                                 0,             size,
                                        0 ); 

    PLA_Obj_horz_split_2( Ainv_BL, size, &Ainv_10_11,
                                          &Ainv_BL );

    PLA_Mscalar_create_conf_to( A_11, PLA_ALL_ROWS, PLA_ALL_COLS,
				 &A_11_msc );
    PLA_Copy( A_11, A_11_msc );

    PLA_Obj_set_orientation( Ainv_10_11, PLA_PROJ_ONTO_ROW );
    PLA_Pmvector_create_conf_to( Ainv_BL, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
				  size, &Ainv_10_11_dup );
    PLA_Pmvector_create_conf_to( Ainv_BL, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
				  size, &A_21_dup );

    /*    PLA_Trsm( PLA_SIDE_LEFT,    PLA_LOWER_TRIANGULAR,
	       PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG,
	       one, A_11, Ainv_10_11 );  */

    PLA_Obj_set_to_identity( A_11 );
    
    PLA_Copy( Ainv_10_11, Ainv_10_11_dup );

    PLA_Local_trsm( PLA_SIDE_LEFT,    PLA_LOWER_TRIANGULAR,
	       PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG,
	       one, A_11_msc, Ainv_10_11_dup ); 

    PLA_Copy( A_21,       A_21_dup );

    PLA_Obj_set_to_zero( A_21 );
 
    PLA_Local_gemm( PLA_NO_TRANS, PLA_NO_TRANS, 
	       minus_one, A_21_dup, Ainv_10_11_dup,
	       one, Ainv_BL );

    PLA_Copy( Ainv_10_11_dup, Ainv_10_11 );
  }

  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_21 );
  PLA_Obj_free( &Ainv );
  PLA_Obj_free( &Ainv_BL );
  PLA_Obj_free( &Ainv_10_11 );
  PLA_Obj_free( &minus_one );
  PLA_Obj_free( &zero );
  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}



