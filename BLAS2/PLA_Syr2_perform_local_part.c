/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Syr2_perform_local_part( int uplo, 
				  PLA_Obj alpha, PLA_Obj xdpmv, 
				                  PLA_Obj ydpmv_trans,
				  PLA_Obj ydpmv, PLA_Obj xdpmv_trans,
                             PLA_Obj A )
{
  int 
    length, local_size, length_1,
    count, me;

  PLA_Obj 
    A11           = NULL,     
    A21           = NULL,       A22           = NULL,
    xdpmv_1       = NULL,       xdpmv_2       = NULL,
    xdpmv_trans_1       = NULL, xdpmv_trans_2       = NULL,
    ydpmv_1       = NULL,       ydpmv_2       = NULL,
    ydpmv_trans_1       = NULL, ydpmv_trans_2       = NULL;


  PLA_Obj_local_width( A, &local_size );
  if ( 0 == local_size ) return;
  PLA_Obj_local_length( A, &local_size );
  if ( 0 == local_size ) return;

  PLA_Obj_global_length( A, &length );

  if ( length >= 256 ) {
    length_1 = (length/2);
    length_1 = length_1 - length_1 % 64 + 64; 
    if ( length_1 >= length ) length = length/2;
    PLA_Obj_split_4 ( A, length_1, length_1,      &A11, PLA_DUMMY,
                                                  &A21, &A22 );
    PLA_Obj_horz_split_2( xdpmv, length_1,       &xdpmv_1,
                                                  &xdpmv_2 );
    PLA_Obj_vert_split_2( ydpmv_trans, length_1,  
                                  &ydpmv_trans_1, &ydpmv_trans_2 );
    PLA_Obj_horz_split_2( ydpmv, length_1,       &ydpmv_1,
                                                  &ydpmv_2 );
    PLA_Obj_vert_split_2( xdpmv_trans, length_1,  
                                  &xdpmv_trans_1, &xdpmv_trans_2 );

    PLA_Local_ger( alpha, xdpmv_2, ydpmv_trans_1, A21 );
    PLA_Local_ger( alpha, ydpmv_2, xdpmv_trans_1, A21 );

    PLA_Syr2_perform_local_part( uplo, alpha, 
				  xdpmv_1, ydpmv_trans_1, 
				  ydpmv_1, xdpmv_trans_1, 
				  A11 );
    PLA_Syr2_perform_local_part( uplo, alpha, 
				  xdpmv_2, ydpmv_trans_2, 
				  ydpmv_2, xdpmv_trans_2, 
				  A22 );
  } 
  else {    
    PLA_Syr2_perform_local_part_by_panels( uplo, alpha, 
					    xdpmv, ydpmv_trans, 
					    ydpmv, xdpmv_trans, 
					    A );
  }

  PLA_Obj_free( &A11 );
  PLA_Obj_free( &A21 );                PLA_Obj_free( &A22 );
  PLA_Obj_free( &xdpmv_1 );            PLA_Obj_free( &xdpmv_2 );
  PLA_Obj_free( &ydpmv_1 );            PLA_Obj_free( &ydpmv_2 );
  PLA_Obj_free( &xdpmv_trans_1 );      PLA_Obj_free( &xdpmv_trans_2 );
  PLA_Obj_free( &ydpmv_trans_1 );      PLA_Obj_free( &ydpmv_trans_2 );

  return( PLA_SUCCESS );
}


int PLA_Syr2_perform_local_part_by_panels( 
                                 int uplo, 
                                 PLA_Obj alpha, 
				 PLA_Obj xdpmv, PLA_Obj ydpmv_trans,
				 PLA_Obj ydpmv, PLA_Obj xdpmv_trans,
                                 PLA_Obj A )
{
  int 
    size, size_left, size_top, dummy, owner_left, mycol, local_size;

  PLA_Obj 
    ABR    = NULL,
    A11    = NULL,   A21     = NULL,
    xdpmv_cur = NULL, xdpmv_top = NULL,
    ydpmv_trans_cur = NULL, ydpmv_trans_left = NULL,
    ydpmv_cur = NULL, ydpmv_top = NULL,
    xdpmv_trans_cur = NULL, xdpmv_trans_left = NULL;

  PLA_Template 
    templ = NULL;

  PLA_Obj_local_width( A, &local_size );
  if ( local_size == 0 ) return;
  PLA_Obj_local_length( A, &local_size );
  if ( local_size == 0 ) return;

  PLA_Obj_view_all( A, &ABR );
  PLA_Obj_view_all( xdpmv, &xdpmv_cur );
  PLA_Obj_view_all( ydpmv, &ydpmv_cur );
  PLA_Obj_view_all( xdpmv_trans, &xdpmv_trans_cur );
  PLA_Obj_view_all( ydpmv_trans, &ydpmv_trans_cur );

  PLA_Obj_template( A, &templ );
  PLA_Temp_comm_row_rank( templ, &mycol );

  while ( TRUE ) {
    PLA_Obj_split_size ( ABR, PLA_SIDE_TOP,  &size_top,  &dummy );
    PLA_Obj_split_size ( ABR, PLA_SIDE_LEFT, &size_left, &owner_left );

    if ( ( size = min ( size_top, size_left ) ) == 0 ) break;

    PLA_Obj_split_4( ABR, size, size,     &A11, PLA_DUMMY,
		                           &A21, &ABR );

    PLA_Obj_horz_split_2( xdpmv_cur, size,  &xdpmv_top,
                                             &xdpmv_cur );
    PLA_Obj_vert_split_2( ydpmv_trans_cur, size, 
                                 &ydpmv_trans_left, &ydpmv_trans_cur );

    PLA_Obj_horz_split_2( ydpmv_cur, size,  &ydpmv_top,
                                             &ydpmv_cur );
    PLA_Obj_vert_split_2( xdpmv_trans_cur, size, 
                                 &xdpmv_trans_left, &xdpmv_trans_cur );

    if ( owner_left == mycol ) {
      PLA_Local_syr2 ( uplo, alpha, xdpmv_top, ydpmv_trans_left, A11 );
      PLA_Local_ger (alpha, xdpmv_cur, ydpmv_trans_left, A21 );
      PLA_Local_ger (alpha, ydpmv_cur, xdpmv_trans_left, A21 );
    }
  }

  PLA_Obj_free( &ABR );
  PLA_Obj_free( &A11 );
  PLA_Obj_free( &A21 );
  PLA_Obj_free( &xdpmv_cur );
  PLA_Obj_free( &xdpmv_top );
  PLA_Obj_free( &ydpmv_trans_cur );
  PLA_Obj_free( &ydpmv_trans_left );
  PLA_Obj_free( &ydpmv_cur );
  PLA_Obj_free( &ydpmv_top );
  PLA_Obj_free( &xdpmv_trans_cur );
  PLA_Obj_free( &xdpmv_trans_left );

  return( PLA_SUCCESS );
}






