/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int nb_out_input = 0;

int PLA_JVTrsv_sub_ln( int, PLA_Obj, PLA_Obj, PLA_Obj, PLA_Obj );

int PLA_JVTrsv_ln( int diag, PLA_Obj A, PLA_Obj b )
{
  int 
    nb_distr, nb_out, nprocs, size;

  PLA_Obj   
    A_BL = NULL,       A_BR = NULL,
    A_10 = NULL,       A_11 = NULL,
    A_20 = NULL,       A_21 = NULL,
    bc = NULL,         xr = NULL,
    b_1 = NULL,        b_B = NULL,
    bc_1 = NULL,       bc_B = NULL,
    xr_R = NULL,
    xr_0 = NULL,       xr_1 = NULL,
    minus_one = NULL,  one = NULL;

  PLA_Template 
    templ = NULL;

  PLA_Obj_template( A, &templ );
  PLA_Temp_nb( templ, &nb_distr );
  PLA_Temp_comm_all_size( templ, &nprocs );
  nb_out = nprocs * nb_distr;   

  PLA_Pmvector_create_conf_to(A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1, &bc );
  PLA_Pmvector_create_conf_to(A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1, &xr );

  PLA_Obj_set_to_zero( bc );
  PLA_Obj_set_to_zero( xr );

  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  PLA_Obj_vert_split_2( A, 0,       &A_BL, &A_BR ); 
  PLA_Obj_vert_split_2( xr, 0,      &xr_0, &xr_R );
  PLA_Obj_view_all( bc, &bc_B );
  PLA_Obj_view_all( b,  &b_B );

  while( TRUE ) {
    PLA_Obj_global_length( A_BR, &size );
    if ( 0 == ( size = min( size, nb_out ) ) ) break;

    PLA_Obj_horz_split_2( A_BL, size,    &A_10,
                                         &A_20 );
    PLA_Obj_split_4( A_BR, size, size,   &A_11, PLA_DUMMY, 
                                         &A_21, &A_BR );    

    PLA_Obj_horz_split_2( b_B, size,     &b_1,
                                         &b_B );
    PLA_Obj_horz_split_2( bc_B, size,    &bc_1,
                                         &bc_B );
    PLA_Obj_vert_split_2( xr_R, size,    &xr_1, &xr_R );

    PLA_Local_gemv( PLA_NO_TRANS, minus_one, A_10, xr_0, one, bc_1 );
    
    PLA_JVTrsv_sub_ln( diag, A_11, b_1, bc_1, xr_1 );

    PLA_Obj_view_shift( A_BL,       size,
                              0,             size,
                                      0 );
    PLA_Obj_view_shift( xr_0,       0,
                              0,             size,
                                    0 );
  }

  PLA_Obj_free(&A_BL);       PLA_Obj_free(&A_BR);
  PLA_Obj_free(&A_10);       PLA_Obj_free(&A_11);
  PLA_Obj_free(&A_20);       PLA_Obj_free(&A_21);

  PLA_Obj_free(&b_1);        PLA_Obj_free(&b_B);
  PLA_Obj_free(&bc);
  PLA_Obj_free(&bc_1);       PLA_Obj_free(&bc_B);
  PLA_Obj_free(&xr);         PLA_Obj_free(&xr_R);   
  PLA_Obj_free(&xr_0);       PLA_Obj_free(&xr_1);

  PLA_Obj_free(&minus_one);  PLA_Obj_free(&one);


  return PLA_SUCCESS;
}



int PLA_JVTrsv_sub_ln( int diag, PLA_Obj A, PLA_Obj b, 
                                   PLA_Obj bc, PLA_Obj xr )
{
  int 
    mycol, myrow, 
    size, size_row, size_col, 
    owner_col, owner_row;

  PLA_Obj   
    A_BL = NULL,       A_BR = NULL,
    A_10 = NULL,       A_11 = NULL,
    A_20 = NULL,       A_21 = NULL,
    b_1 = NULL,        b_B = NULL,
    bc_1 = NULL,       bc_B = NULL,
    xr_R = NULL,
    xr_0 = NULL,       xr_1 = NULL,
    minus_one = NULL,  one = NULL;

  PLA_Template 
    templ = NULL;

  PLA_Obj_template( A, &templ );
  PLA_Temp_comm_row_rank( templ, &mycol );
  PLA_Temp_comm_col_rank( templ, &myrow );

  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  PLA_Obj_vert_split_2( A, 0,       &A_BL, &A_BR ); 
  PLA_Obj_vert_split_2( xr, 0,      &xr_0, &xr_R );
  PLA_Obj_view_all( bc, &bc_B );
  PLA_Obj_view_all( b,  &b_B );

  while( TRUE ) {
    PLA_Obj_split_size( A_BR, PLA_SIDE_TOP, &size_row, &owner_row );
    PLA_Obj_split_size( A_BR, PLA_SIDE_LEFT, &size_col, &owner_col );
    if ( 0 == ( size = min( size_row, size_col ) ) ) break;

    PLA_Obj_horz_split_2( A_BL, size,    &A_10,
                                         &A_20 );
    PLA_Obj_split_4( A_BR, size, size, &A_11, PLA_DUMMY, 
                                       &A_21, &A_BR );    

    PLA_Obj_horz_split_2( b_B, size,     &b_1,
                                         &b_B );
    PLA_Obj_horz_split_2( bc_B, size,    &bc_1,
                                         &bc_B );
    PLA_Obj_vert_split_2( xr_R, size,    &xr_1, &xr_R );

    {
      int local_width;
      PLA_Obj_local_width( xr_0, &local_width );
      if ( myrow == owner_row && local_width > 0 ) 
	PLA_Local_gemv( PLA_NO_TRANS, minus_one, A_10, xr_0, one, bc_1 );
    }

/*    if (myrow == owner_row ) */
/*      PLA_Reduce_x( PLA_SHAPE_GENERAL, bc_1, one, b_1 ); */
    PLA_Reduce_add( bc_1, MPI_SUM, b_1 );

    if (myrow == owner_row && mycol == owner_col ) 
      PLA_Local_trsv( PLA_LOWER_TRIANGULAR, PLA_NO_TRANS, diag,
		      A_11, b_1 );
    
/*    if ( mycol == owner_col ) */
      PLA_Copy( b_1, xr_1 );

    PLA_Obj_view_shift( A_BL,       size,
                              0,             size,
                                      0 );
    PLA_Obj_view_shift( xr_0,       0,
                              0,             size,
                                    0 );
  }
  PLA_Obj_free(&A_BL);       PLA_Obj_free(&A_BR);
  PLA_Obj_free(&A_10);       PLA_Obj_free(&A_11);
  PLA_Obj_free(&A_20);       PLA_Obj_free(&A_21);

  PLA_Obj_free(&b_1);        PLA_Obj_free(&b_B);
  PLA_Obj_free(&bc_1);       PLA_Obj_free(&bc_B);
  PLA_Obj_free(&xr_R);   
  PLA_Obj_free(&xr_0);       PLA_Obj_free(&xr_1);

  PLA_Obj_free(&minus_one);  PLA_Obj_free(&one);

  return PLA_SUCCESS;
}

