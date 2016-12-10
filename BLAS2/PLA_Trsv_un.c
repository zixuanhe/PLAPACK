/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

void PLA_JVTrsv_sub_un( int, PLA_Obj, PLA_Obj, PLA_Obj, PLA_Obj );

void PLA_JVTrsv_un( int diag, PLA_Obj A, PLA_Obj b )
{
  PLA_Obj   A_TL = NULL,       A_TR = NULL,
            A_01 = NULL,       A_02 = NULL,
            A_11 = NULL,       A_12 = NULL,
            bc = NULL,         xr = NULL,
            b_1 = NULL,        b_T = NULL,
            bc_1 = NULL,       bc_T = NULL,
            xr_L = NULL,
            xr_1 = NULL,       xr_2 = NULL,
            minus_one = NULL,  one = NULL;
  PLA_Template templ = NULL;
  int nb_distr, nb_out, nprocs, size;

  PLA_Obj_template( A, &templ );
  PLA_Temp_nb( templ, &nb_distr );
  PLA_Temp_comm_all_size( templ, &nprocs );
  nb_out = nprocs * nb_distr;   

  PLA_Pmvector_create_conf_to(A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1, &bc );
  PLA_Pmvector_create_conf_to(A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1, &xr );

  PLA_Obj_set_to_zero( bc );
  PLA_Obj_set_to_zero( xr );

  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  PLA_Obj_global_length( A, &size );
  PLA_Obj_vert_split_2( A, size,       &A_TL, &A_TR ); 
  PLA_Obj_vert_split_2( xr, size,      &xr_L, &xr_2 );
  PLA_Obj_view_all( bc, &bc_T );
  PLA_Obj_view_all( b,  &b_T );

  while( TRUE ) {
    PLA_Obj_global_length( A_TL, &size );
    if ( 0 == ( size = min( size, nb_out ) ) ) break;

    PLA_Obj_horz_split_2( A_TR, -size,    &A_02,
                                          &A_12 );
    PLA_Obj_split_4( A_TL, -size, -size, &A_TL,     &A_01, 
                                         PLA_DUMMY, &A_11 );    

    PLA_Obj_horz_split_2( b_T, -size,     &b_T,
                                          &b_1 );
    PLA_Obj_horz_split_2( bc_T, -size,    &bc_T,
                                          &bc_1 );
    PLA_Obj_vert_split_2( xr_L, -size,    &xr_L, &xr_1 );

    PLA_Local_gemv( PLA_NO_TRANS, minus_one, A_12, xr_2, one, bc_1 );

    PLA_JVTrsv_sub_un( diag, A_11, b_1, bc_1, xr_1 );

    PLA_Obj_view_shift( A_TR,        0,
                              -size,       0,
                                    -size );
    PLA_Obj_view_shift( xr_2,       0,
                              -size,     0,
                                    0 );
  }

  PLA_Obj_free(&A_TL);       PLA_Obj_free(&A_TR);
  PLA_Obj_free(&A_01);       PLA_Obj_free(&A_02);
  PLA_Obj_free(&A_11);       PLA_Obj_free(&A_12);

  PLA_Obj_free(&b_1);        PLA_Obj_free(&b_T);
  PLA_Obj_free(&bc);
  PLA_Obj_free(&bc_1);       PLA_Obj_free(&bc_T);
  PLA_Obj_free(&xr);         PLA_Obj_free(&xr_L);   
  PLA_Obj_free(&xr_1);       PLA_Obj_free(&xr_2);

  PLA_Obj_free(&minus_one);  PLA_Obj_free(&one);


  return;
}



void PLA_JVTrsv_sub_un( int diag, PLA_Obj A, PLA_Obj b, 
                                   PLA_Obj bc, PLA_Obj xr )
{
  PLA_Obj   A_TL = NULL,       A_TR = NULL,
            A_01 = NULL,       A_02 = NULL,
            A_11 = NULL,       A_12 = NULL,
            b_1 = NULL,        b_T = NULL,
            bc_1 = NULL,       bc_T = NULL,
            xr_L = NULL,
            xr_1 = NULL,       xr_2 = NULL,
            minus_one = NULL,  one = NULL;
  PLA_Template templ = NULL;
  int mycol, myrow, size, size_row, size_col, owner_col, owner_row;

  PLA_Obj_template( A, &templ );
  PLA_Temp_comm_row_rank( templ, &mycol );
  PLA_Temp_comm_col_rank( templ, &myrow );

  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  PLA_Obj_global_length( A, &size );
  PLA_Obj_vert_split_2( A, size,       &A_TL, &A_TR ); 
  PLA_Obj_vert_split_2( xr, size,      &xr_L, &xr_2 );
  PLA_Obj_view_all( bc, &bc_T );
  PLA_Obj_view_all( b,  &b_T );

  while( TRUE ) {
    PLA_Obj_split_size( A_TL, PLA_SIDE_BOTTOM, &size_row, &owner_row );
    PLA_Obj_split_size( A_TL, PLA_SIDE_RIGHT,  &size_col, &owner_col );
    if ( 0 == ( size = min( size_row, size_col ) ) ) break;

    PLA_Obj_horz_split_2( A_TR, -size,    &A_02,
                                          &A_12 );
    PLA_Obj_split_4( A_TL, -size, -size, &A_TL,     &A_01,
                                         PLA_DUMMY, &A_11 );    

    PLA_Obj_horz_split_2( b_T, -size,     &b_T,
                                          &b_1 );
    PLA_Obj_horz_split_2( bc_T, -size,    &bc_T,
                                          &bc_1 );
    PLA_Obj_vert_split_2( xr_L, -size,    &xr_L, &xr_1 );

    {
      int local_width;
      PLA_Obj_local_width( xr_2, &local_width );
      if ( myrow == owner_row && local_width > 0 ) 
	PLA_Local_gemv( PLA_NO_TRANS, minus_one, A_12, xr_2, one, bc_1 );
    }

/*    if (myrow == owner_row ) */
/*      PLA_Reduce_x( PLA_SHAPE_GENERAL, bc_1, one, b_1 ); */
    PLA_Reduce_add( bc_1, MPI_SUM, b_1 );

    if (myrow == owner_row && mycol == owner_col ) 
      PLA_Local_trsv( PLA_UPPER_TRIANGULAR, PLA_NO_TRANS, diag,
		      A_11, b_1 );
    
/*     if ( mycol == owner_col ) */
      PLA_Copy( b_1, xr_1 );

    PLA_Obj_view_shift( A_TR,         0,
                              -size,     0,
                                      -size );
    PLA_Obj_view_shift( xr_2,       0,
                              -size,      0,
                                    0 );
  }
  PLA_Obj_free(&A_TL);       PLA_Obj_free(&A_TR);
  PLA_Obj_free(&A_01);       PLA_Obj_free(&A_02);
  PLA_Obj_free(&A_11);       PLA_Obj_free(&A_12);

  PLA_Obj_free(&b_1);        PLA_Obj_free(&b_T);
  PLA_Obj_free(&bc_1);       PLA_Obj_free(&bc_T);
  PLA_Obj_free(&xr_L);   
  PLA_Obj_free(&xr_1);       PLA_Obj_free(&xr_2);

  PLA_Obj_free(&minus_one);  PLA_Obj_free(&one);

  return;
}

