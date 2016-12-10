/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Symmetric_invert_lower( PLA_Obj A )
{
  int       
    value = 0,
    size, k, nb_alg;   

  PLA_Obj  
    ABR = NULL, A11 = NULL, A11_dup = NULL,
    ATL = NULL,
    ABL = NULL, A10 = NULL, A21 = NULL,
    A1_mv = NULL, A1_mv_T = NULL, A1_mv_1 = NULL, A1_mv_B = NULL,
    A1_dpmv = NULL, A1_dpmv_T = NULL, A1_dpmv_1 = NULL, A1_dpmv_B = NULL,
    A1_T_dpmv = NULL, A1_T_dpmv_L = NULL, A1_T_dpmv_1 = NULL, 
    A1_T_dpmv_R = NULL,
    minus_one = NULL, zero = NULL, one = NULL;

  PLA_Template
    templ;

  PLA_Obj_template( A, &templ );
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  /* Create usual duplicated scalar constants */
  PLA_Create_constants_conf_to( A, &minus_one, &zero, &one );

  /* View A = / ATL || ATR \ 
              | ====  ==== | 	      where ATL is 0x0 
              \ ABL || ABR /                                 */

  PLA_Obj_split_4( A, 0, 0, &ATL, PLA_DUMMY,
		             &ABL, &ABR         );

  PLA_Mvector_create_conf_to( A, nb_alg, &A1_mv );
  PLA_Obj_horz_split_2( A1_mv, 0, &A1_mv_T,
			           &A1_mv_B );

  PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
				nb_alg, &A1_dpmv );
  PLA_Obj_horz_split_2( A1_dpmv, 0, &A1_dpmv_T,
			           &A1_dpmv_B );

  PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
				nb_alg, &A1_T_dpmv );
  PLA_Obj_vert_split_2( A1_T_dpmv, 0, &A1_T_dpmv_L, &A1_T_dpmv_R );

  k = 0;
  while ( TRUE ) {
    /* Determine size of current panel */
    PLA_Obj_global_length( ABR, &size );
    if ( 0 == ( size = min( size, nb_alg ) ) ) break;

    PLA_Obj_horz_split_2( ABL, size,       &A10,
                                            PLA_DUMMY );
    PLA_Obj_set_orientation( A10, PLA_PROJ_ONTO_ROW );

    PLA_Obj_split_4(      ABR, size, size,        &A11, PLA_DUMMY,
                                                   &A21, &ABR );

    PLA_Mscalar_create_conf_to( A11, PLA_ALL_ROWS, PLA_ALL_COLS, 
				 &A11_dup );

    PLA_Copy( A11, A11_dup );

    value = PLA_Local_chol( PLA_LOWER_TRIANGULAR, A11_dup );
    if ( value != 0 ) {
      value = k + value;
      break;
    }

    if ( size < nb_alg ){
      PLA_Obj_vert_split_2( A1_mv, size, &A1_mv, PLA_DUMMY );
      PLA_Obj_vert_split_2( A1_mv_T, size, &A1_mv_T, PLA_DUMMY );
      PLA_Obj_vert_split_2( A1_mv_B, size, &A1_mv_B, PLA_DUMMY );

      PLA_Obj_vert_split_2( A1_dpmv, size,   &A1_dpmv, PLA_DUMMY );
      PLA_Obj_vert_split_2( A1_dpmv_T, size, &A1_dpmv_T, PLA_DUMMY );
      PLA_Obj_vert_split_2( A1_dpmv_B, size, &A1_dpmv_B, PLA_DUMMY );

      PLA_Obj_horz_split_2( A1_T_dpmv, size,   &A1_T_dpmv, 
			                        PLA_DUMMY );
      PLA_Obj_horz_split_2( A1_T_dpmv_L, size, &A1_T_dpmv_L, 
			                        PLA_DUMMY );
      PLA_Obj_horz_split_2( A1_T_dpmv_R, size, &A1_T_dpmv_R, 
			                        PLA_DUMMY );
    }

    PLA_Obj_horz_split_2( A1_mv_B, size, &A1_mv_1,
                                          &A1_mv_B );

    PLA_Obj_horz_split_2( A1_dpmv_B, size, &A1_dpmv_1,
                                            &A1_dpmv_B );

    PLA_Obj_vert_split_2( A1_T_dpmv_R, size, &A1_T_dpmv_1, &A1_T_dpmv_R );

    PLA_Copy( A10, A1_mv_T );
    PLA_Copy( A21, A1_mv_B );
    
    PLA_Obj_set_to_identity( A1_mv_1 );


    PLA_Local_trsm( PLA_SIDE_RIGHT, PLA_LOWER_TRIANGULAR,
		     PLA_TRANSPOSE,  PLA_NONUNIT_DIAG,
		     one, A11_dup, A1_mv );

    PLA_Copy( A1_mv, A1_dpmv );
    PLA_Copy( A1_mv, A1_T_dpmv );

    PLA_Obj_set_to_zero( A10 );
    PLA_Obj_set_to_zero( A11 );
    PLA_Obj_set_to_zero( A21 );

    PLA_Obj_view_shift( A1_mv_T,    0,
                              0,           0,
                                    size );

    PLA_Obj_view_shift( ABL,        size,
                              0,           size,
                                       0  );

    PLA_Obj_view_shift( ATL,          0,
                              0,           size,
                                       size  );

    PLA_Obj_view_shift( A1_dpmv_T,    0,
                              0,           0,
                                       size  );

    PLA_Obj_view_shift( A1_T_dpmv_L,    0,
                              0,             size,
                                         0  );

    PLA_Syrk_perform_local_part( PLA_LOWER_TRIANGULAR, 
				  one, A1_dpmv_T, A1_T_dpmv_L,
				  one, ATL );

    PLA_Local_gemm( PLA_NO_TRANS, PLA_NO_TRANS, 
		     minus_one, A1_dpmv_B, A1_T_dpmv_L, one, ABL );

    PLA_Syrk_perform_local_part( PLA_LOWER_TRIANGULAR, 
				  minus_one, A1_dpmv_B, A1_T_dpmv_R,
				  one, ABR );

    k = k+size;
  }

  PLA_Obj_free( &ABR );
  PLA_Obj_free( &A11 );
  PLA_Obj_free( &A11_dup );
  PLA_Obj_free( &ATL );
  PLA_Obj_free( &ABL );
  PLA_Obj_free( &A10 );
  PLA_Obj_free( &A21 );
  PLA_Obj_free( &A1_mv );
  PLA_Obj_free( &A1_mv_T );
  PLA_Obj_free( &A1_mv_1 );
  PLA_Obj_free( &A1_mv_B );
  PLA_Obj_free( &A1_dpmv );
  PLA_Obj_free( &A1_dpmv_T );
  PLA_Obj_free( &A1_dpmv_1 );
  PLA_Obj_free( &A1_dpmv_B );
  PLA_Obj_free( &A1_T_dpmv );
  PLA_Obj_free( &A1_T_dpmv_L );
  PLA_Obj_free( &A1_T_dpmv_1 );
  PLA_Obj_free( &A1_T_dpmv_R );
  PLA_Obj_free( &minus_one );
  PLA_Obj_free( &zero );
  PLA_Obj_free( &one );

  return value;
}
