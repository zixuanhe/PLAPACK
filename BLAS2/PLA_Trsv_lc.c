/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

void PLA_JVTrsv_sub_lc( int, PLA_Obj, PLA_Obj, PLA_Obj, PLA_Obj );

void PLA_JVTrsv_lc( int diag, PLA_Obj A, PLA_Obj b )
{
  PLA_Obj   A_TL = NULL,       A_BL = NULL,
            A_10 = NULL,       A_11 = NULL,
            A_20 = NULL,       A_21 = NULL,
            br = NULL,         xc = NULL,
            b_1 = NULL,        b_L = NULL,
            br_1 = NULL,       br_L = NULL,
            xc_T = NULL,
            xc_1 = NULL,       xc_2 = NULL,
            minus_one = NULL,  one = NULL;
  PLA_Template templ = NULL;
  int nb_distr, nb_out, nprocs, size;

  PLA_Obj_template( A, &templ );
  PLA_Temp_nb( templ, &nb_distr );
  PLA_Temp_comm_all_size( templ, &nprocs );
  nb_out = nprocs * nb_distr;   
/*  nb_out = 10000; */

  PLA_Pmvector_create_conf_to(A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1, &xc );
  PLA_Pmvector_create_conf_to(A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1, &br );

  PLA_Obj_set_to_zero( br );
  PLA_Obj_set_to_zero( xc );

  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  PLA_Obj_global_length( A, &size );
  PLA_Obj_horz_split_2( A, size,       &A_TL, 
                                       &A_BL ); 
  PLA_Obj_horz_split_2( xc, size,      &xc_T, 
                                       &xc_2 );
  PLA_Obj_view_all( br, &br_L );
  PLA_Obj_view_all( b,  &b_L );

  while( TRUE ) {
    PLA_Obj_global_length( A_TL, &size );
    if ( 0 == ( size = min( size, nb_out ) ) ) break;

    PLA_Obj_vert_split_2( A_BL, -size,   &A_20, &A_21 );
    PLA_Obj_split_4( A_TL, -size, -size,   &A_TL, PLA_DUMMY, 
                                           &A_10, &A_11 );    

    PLA_Obj_horz_split_2( b_L, -size,    &b_L,  
                                         &b_1 );
    PLA_Obj_vert_split_2( br_L, -size,   &br_L, &br_1 );
    PLA_Obj_horz_split_2( xc_T, -size,   &xc_T, &xc_1 );

    PLA_Local_gemv( PLA_CONJ_TRANS, minus_one, A_21, xc_2, one, br_1 ); 

    PLA_JVTrsv_sub_lc( diag, A_11, b_1, br_1, xc_1 ); 

    PLA_Obj_view_shift( A_BL,       -size,
                              0,            -size,
                                       0 );
    PLA_Obj_view_shift( xc_2,      -size,
                              0,             0,
                                    0 );
  }

  PLA_Obj_free(&A_TL);       PLA_Obj_free(&A_BL);
  PLA_Obj_free(&A_10);       PLA_Obj_free(&A_11);
  PLA_Obj_free(&A_20);       PLA_Obj_free(&A_21);

  PLA_Obj_free(&b_1);        PLA_Obj_free(&b_L);
  PLA_Obj_free(&br);
  PLA_Obj_free(&br_1);       PLA_Obj_free(&br_L);
  PLA_Obj_free(&xc);         PLA_Obj_free(&xc_T);   
  PLA_Obj_free(&xc_1);       PLA_Obj_free(&xc_2);

  PLA_Obj_free(&minus_one);  PLA_Obj_free(&one);


  return;
}



void PLA_JVTrsv_sub_lc( int diag, PLA_Obj A, PLA_Obj b, 
                                   PLA_Obj br, PLA_Obj xc )
{
  PLA_Obj   A_TL = NULL,       A_BL = NULL,
            A_10 = NULL,       A_11 = NULL,
            A_20 = NULL,       A_21 = NULL,
            b_1 = NULL,        b_L = NULL,
            br_1 = NULL,       br_L = NULL,
            xc_T = NULL,
            xc_1 = NULL,       xc_2 = NULL,
            minus_one = NULL,  one = NULL;
  PLA_Template templ = NULL;
  int mycol, myrow, size, size_row, size_col, owner_col, owner_row;

  PLA_Obj_template( A, &templ );
  PLA_Temp_comm_row_rank( templ, &mycol );
  PLA_Temp_comm_col_rank( templ, &myrow );

  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  PLA_Obj_global_length( A, &size );
  PLA_Obj_horz_split_2( A, size,       &A_TL, 
                                       &A_BL ); 
  PLA_Obj_horz_split_2( xc, size,      &xc_T, 
                                       &xc_2 );
  PLA_Obj_view_all( br, &br_L );
  PLA_Obj_view_all( b,  &b_L );
			     
  while( TRUE ) {
    PLA_Obj_split_size( A_TL, PLA_SIDE_BOTTOM, &size_row, &owner_row );
    PLA_Obj_split_size( A_TL, PLA_SIDE_RIGHT,  &size_col, &owner_col );
    if ( 0 == ( size = min( size_row, size_col ) ) ) break;

    PLA_Obj_vert_split_2( A_BL, -size,   &A_20, &A_21 );
    PLA_Obj_split_4( A_TL, -size, -size, &A_TL, PLA_DUMMY, 
                                         &A_10, &A_11 );    

    PLA_Obj_horz_split_2( b_L, -size,    &b_L,
                                         &b_1 );
    PLA_Obj_vert_split_2( br_L, -size,   &br_L, &br_1 );
    PLA_Obj_horz_split_2( xc_T, -size,   &xc_T,  
                                         &xc_1 );

    {
      int local_length;
      PLA_Obj_local_length( xc_2, &local_length );
      if ( mycol == owner_col && local_length > 0 ) 
	PLA_Local_gemv( PLA_CONJ_TRANS, minus_one, A_21, xc_2, one, br_1 );
    }

/*    if (mycol == owner_col ) */{
/*      PLA_Reduce_x( PLA_SHAPE_GENERAL, br_1, one, b_1 ); */
      PLA_Reduce_add( br_1, MPI_SUM, b_1 );
    }
    if (mycol == owner_col && myrow == owner_row ) 
      PLA_Local_trsv( PLA_LOWER_TRIANGULAR, PLA_CONJ_TRANS, diag,
		      A_11, b_1 );
    
/*    if ( myrow == owner_row ) */{
      PLA_Copy( b_1, xc_1 );
    }

    PLA_Obj_view_shift( A_BL,         -size,
                              0,             -size,
                                        0 );
    PLA_Obj_view_shift( xc_2,       -size,
                              0,             0,
                                    0 );
  }
  PLA_Obj_free(&A_TL);       PLA_Obj_free(&A_BL);
  PLA_Obj_free(&A_10);       PLA_Obj_free(&A_11);
  PLA_Obj_free(&A_20);       PLA_Obj_free(&A_21);

  PLA_Obj_free(&b_1);        PLA_Obj_free(&b_L);
  PLA_Obj_free(&br_1);       PLA_Obj_free(&br_L);
  PLA_Obj_free(&xc_T);   
  PLA_Obj_free(&xc_1);       PLA_Obj_free(&xc_2);

  PLA_Obj_free(&minus_one);  PLA_Obj_free(&one);

  return;
}


