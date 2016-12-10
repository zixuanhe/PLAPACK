/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Triangular_Lt_L( PLA_Obj A )
{
  int 
    nb_alg, size, value = PLA_SUCCESS;

  PLA_Obj
    A_TL = NULL, A_BL = NULL, A_10_11 = NULL, A_11 = NULL, 
    A_10_11_dpmv_rows = NULL, A_10_11_dpmv_cols = NULL,
    one = NULL;

  PLA_Template
    templ;

  PLA_Obj_template( A, &templ );
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  PLA_Create_constants_conf_to( A, NULL, NULL, &one );

  PLA_Obj_split_4( A, 0, 0, &A_TL, PLA_DUMMY,
		             &A_BL, PLA_DUMMY );

  while( TRUE ){
    PLA_Obj_global_length( A_BL, &size );
    if ( 0 == ( size = min( size, nb_alg ) ) ) break;

    PLA_Obj_view_shift( A_TL,         0,
                               0,            size,
                                     size );

    PLA_Obj_view_shift( A_BL,         0,
                               0,            size,
                                       0 );

    PLA_Obj_horz_split_2( A_BL, size,   &A_10_11,
			                 &A_BL );

    PLA_Obj_vert_split_2( A_10_11, -size, PLA_DUMMY, &A_11 );

    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, A_11 );

    PLA_Obj_set_orientation( A_10_11, PLA_PROJ_ONTO_ROW );

    PLA_Pmvector_create_conf_to( A_TL,  PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
				  size, &A_10_11_dpmv_cols );

    PLA_Pmvector_create_conf_to( A_TL,  PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
				  size, &A_10_11_dpmv_rows );

    PLA_Copy( A_10_11,           A_10_11_dpmv_rows );

    /*
    {
      PLA_Obj
	A_10_11_mv = NULL;

      PLA_Mvector_create_conf_to( A_10_11, 1, &A_10_11_mv );
      PLA_Copy( A_10_11_dpmv_rows, A_10_11_mv );
      PLA_Copy( A_10_11_mv, A_10_11_dpmv_cols );

      PLA_Obj_free( &A_10_11_mv );
    }
    */
    PLA_Copy( A_10_11,          A_10_11_dpmv_cols );
    PLA_Obj_set_to_zero( A_10_11 );

    PLA_Syrk_perform_local_part( PLA_LOWER_TRIANGULAR, 
				  one, A_10_11_dpmv_cols, A_10_11_dpmv_rows, 
                                  one, A_TL );
  }

  PLA_Obj_free( &A_TL );
  PLA_Obj_free( &A_BL );
  PLA_Obj_free( &A_10_11 );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_10_11_dpmv_rows );
  PLA_Obj_free( &A_10_11_dpmv_cols );
  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}



