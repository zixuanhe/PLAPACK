/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_General_invert_stable( PLA_Obj A )
{
  int
    k, length, align, size, nb_alg, value = PLA_SUCCESS;

  PLA_Obj
    A_R = NULL,       A_1 = NULL,      BA_B = NULL,     BA_1 = NULL,
    A_BR = NULL,      A_11 = NULL,     A_11_msc = NULL,
    A_1_mv = NULL,    A_01_mv = NULL,  A_11_mv = NULL,  A_21_mv = NULL, 
    A_1_mv_B = NULL,  A_1_mv_T = NULL, A_1_dpmv = NULL, 
    BA_1_dpmv = NULL, BA_1_T_mv = NULL,
    pivots = NULL,    pivots_B = NULL, pivots_1 = NULL, pivots_msc = NULL,
    minus_one = NULL, one = NULL;
  PLA_Template
    templ;
    
  k = 0;
  /* Create a vector in which to keep pivot information */
  PLA_Obj_template( A, &templ );
  PLA_Obj_global_length( A, &length );
  PLA_Obj_global_align_col( A, &align );

  PLA_Mvector_create( MPI_INT, length, 1, templ, align, &pivots );

  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  PLA_Obj_view_all( A, &A_R );
  PLA_Obj_view_all( A, &BA_B );
  PLA_Obj_view_all( A, &A_BR );
  PLA_Obj_view_all( pivots, &pivots_B );

  while( TRUE ){
    PLA_Obj_global_width( A_R, &size );
    if ( 0 == ( size = min( size, nb_alg ) ) ) break;

    /* Split off next column panel of updated matrix A */
    PLA_Obj_vert_split_2( A_R, size, &A_1, &A_R );

    /* Copy current column panel to a multivector */
    PLA_Mvector_create_conf_to( A, size, &A_1_mv );
    PLA_Copy( A_1, A_1_mv );

    /* View multivector as / A_01_mv \
                           | A_11_mv |
                           \ A_21_mv / where A_11_mv is the diagonal block 

       Also, let A_1_mv_B = / A_11_mv \
                            \ A_21_mv /, the part to be factored using LU */
    PLA_Obj_horz_split_2( A_1_mv, k,      &A_01_mv,
                                           &A_1_mv_B );
    PLA_Obj_horz_split_2( A_1_mv_B, size, &A_11_mv,
                                           &A_21_mv );

    /* View the next part of the pivot vector */
    PLA_Obj_horz_split_2( pivots_B, size, &pivots_1,
                                           &pivots_B );
    
    /* Create a duplicated multivector in which to store a copy
       of the LU factorization of A_11 (after pivoting) */
    PLA_Mscalar_create_conf_to( A_11_mv,  PLA_ALL_ROWS, PLA_ALL_COLS,
				 &A_11_msc );

    /* Create a duplicated multivector in which to compute the next pivots */
    PLA_Mscalar_create_conf_to( pivots_1, PLA_ALL_ROWS, PLA_ALL_COLS,
				 &pivots_msc );

    /* Factor  / A_11_mv \   using LU factorization with pivoting.
               \ A_21_mv /   A copy of the LU factors of A_11 are
                             returned in A_11_msc */
    value = PLA_LU_right_mv( A_1_mv_B, A_11_msc, pivots_msc );

    if ( value != PLA_SUCCESS ) break;

    /* Store the pivots for future use */
    PLA_Copy( pivots_msc, pivots_1 );

    /* Next we want to update the multivector with
       / - A_01_mv inv(U_11) \
       |   inv( U_11 )       |
       \ - A_21_mv           /
       Notice that A_21_mv has already been updated with 
       L_21_mv = A_21_mv inv( U_11 )
       Step 1: compute / A_01_mv \  inv( U_11 )
                       \    I    /                  */
    PLA_Obj_set_to_identity( A_11_mv );

    PLA_Obj_horz_split_2( A_1_mv, k+size, &A_1_mv_T,
			                   PLA_DUMMY );

    PLA_Local_trsm( PLA_SIDE_RIGHT, PLA_UPPER_TRIANGULAR,
		     PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG,
		     one, A_11_msc, A_1_mv_T );

    /* Step 2: Negate A_01 and A_21                              */
    PLA_Local_scal( minus_one, A_01_mv );
    PLA_Local_scal( minus_one, A_21_mv );

    /* Apply pivots to BA_B */
    PLA_Apply_pivots_to_rows( BA_B, pivots_msc );

    /* Next, update  ( B_10 || A_11 | A_12 ) with 
                     ( inv(L_11) B_10 || inv(L_11) A_11 | inv(L_11) A_12 ) */

   /* Split off next row panel of updated matrix ( B || A ) */
    PLA_Obj_horz_split_2( BA_B, size, &BA_1, 
                                       &BA_B );
    PLA_Obj_set_orientation( BA_1, PLA_PROJ_ONTO_ROW );

    /* Copy current row panel to a multivector.  
       NOTICE: THE COPY TRANSPOSES THE DATA! */
    PLA_Mvector_create_conf_to( A, size, &BA_1_T_mv );
    PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
				  size, &BA_1_dpmv );

    PLA_Obj_split_4(             A_BR, size, size,  &A_11,      PLA_DUMMY,
                                                     PLA_DUMMY, &A_BR );

    PLA_Obj_set_to_identity( A_11 );

    PLA_Copy( BA_1, BA_1_T_mv );

    /* BA_1 <- inv(L_11) BA_1 
       thus
       BA_1_T_mv <- BA_1_T_mv inv(L_11^T) */

    PLA_Local_trsm( PLA_SIDE_RIGHT, PLA_LOWER_TRIANGULAR,
		     PLA_TRANSPOSE, PLA_UNIT_DIAG,
		     one, A_11_msc, BA_1_T_mv );

    /* Copy BA_1_mv and A_1_mv to duplicated multivectors */
    PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
				  size, &A_1_dpmv );
    PLA_Copy( A_1_mv, A_1_dpmv );
    PLA_Copy( BA_1_T_mv, BA_1_dpmv );

    /* Zero out A_1 and BA_1 */

    PLA_Obj_set_to_zero( A_1 );
    PLA_Obj_set_to_zero( BA_1 );

    /*
      Update  local part of
     / B_TL || A_01 | A_02 \    / B_TL || 0 | A_02 \   
     | ====    =========== |    | ====    ======== |   
     | B_10 || A_11 | A_12 | <- |  0   || 0 |   0  | + A_1_mv * BA_1_dpmv
     | ----    ----------- |    | ----    -------- |   
     \ B_20 || A_21 | A_BR /    \ B_20 || 0 | A_BR /   
    */

    PLA_Local_gemm( PLA_NO_TRANS, PLA_NO_TRANS, 
		     one, A_1_dpmv, BA_1_dpmv, one, A );

    k = k+size;
  }

  if ( value == PLA_SUCCESS )
    PLA_Apply_pivots_to_columns_in_reverse( A, pivots );

  PLA_Obj_free( &A_R );               PLA_Obj_free( &A_1 );
  PLA_Obj_free( &BA_B );              PLA_Obj_free( &BA_1 );
  PLA_Obj_free( &A_11 );              PLA_Obj_free( &A_11_msc );
  PLA_Obj_free( &A_11_mv );           PLA_Obj_free( &A_21_mv );
  PLA_Obj_free( &A_1_mv );            PLA_Obj_free( &A_01_mv );
  PLA_Obj_free( &A_1_mv_B );          PLA_Obj_free( &A_1_mv_T );
  PLA_Obj_free( &A_1_dpmv );          PLA_Obj_free( &BA_1_dpmv );         
  PLA_Obj_free( &BA_1_T_mv );         PLA_Obj_free( &A_BR );
  PLA_Obj_free( &pivots );            PLA_Obj_free( &pivots_B );
  PLA_Obj_free( &pivots_1 );          PLA_Obj_free( &pivots_msc );
  PLA_Obj_free( &minus_one );         PLA_Obj_free( &one );

  return value;
}


