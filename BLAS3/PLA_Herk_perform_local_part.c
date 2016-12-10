/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Herk_perform_local_part( int uplo, 
                                  PLA_Obj alpha, PLA_Obj Xdpmv, 
                                                  PLA_Obj Xdpmv_conj_trans,
                                  PLA_Obj beta,  PLA_Obj A )
{
  int 
    length, local_size, length_1,
    count, me;

  PLA_Obj 
    A11           = NULL,     
    A21           = NULL,       A22           = NULL,
    Xdpmv_1       = NULL,       Xdpmv_2       = NULL,
    Xdpmv_conj_trans_1 = NULL,       Xdpmv_conj_trans_2 = NULL;


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
    PLA_Obj_horz_split_2( Xdpmv, length_1,        &Xdpmv_1,
                                                  &Xdpmv_2 );
    PLA_Obj_vert_split_2( Xdpmv_conj_trans, length_1,  
                                  &Xdpmv_conj_trans_1, &Xdpmv_conj_trans_2 );

    PLA_Local_gemm( PLA_NO_TRANS, PLA_NO_TRANS, 
		    alpha, Xdpmv_2, Xdpmv_conj_trans_1, beta, A21 );

    PLA_Herk_perform_local_part( PLA_LOWER_TRIANGULAR, 
                                 alpha, Xdpmv_1, Xdpmv_conj_trans_1,
                                 beta, A11 );
    PLA_Herk_perform_local_part( PLA_LOWER_TRIANGULAR, 
                                 alpha, Xdpmv_2, Xdpmv_conj_trans_2,
                                 beta, A22 );
  } 
  else {    
    PLA_Herk_perform_local_part_by_panels( 
                                 uplo, alpha, Xdpmv, Xdpmv_conj_trans,
                                       beta,  A );
  }
  PLA_Obj_free( &A11 );
  PLA_Obj_free( &A21 );                PLA_Obj_free( &A22 );
  PLA_Obj_free( &Xdpmv_1 );            PLA_Obj_free( &Xdpmv_2 );
  PLA_Obj_free( &Xdpmv_conj_trans_1 );      PLA_Obj_free( &Xdpmv_conj_trans_2 );

  return( PLA_SUCCESS );
}


int PLA_Herk_perform_local_part_by_panels( 
                                 int uplo, 
                                 PLA_Obj alpha, PLA_Obj Xdpmv, 
                                                 PLA_Obj Xdpmv_conj_trans,
                                 PLA_Obj beta,  PLA_Obj A )
{
  int 
    size, size_left, size_top, dummy, owner_left, mycol, local_size;

  PLA_Obj 
    ABR    = NULL,
    A11    = NULL,   A21     = NULL,
    Xdpmv_cur = NULL, Xdpmv_top = NULL,
    Xdpmv_conj_trans_cur = NULL, Xdpmv_conj_trans_left = NULL;

  PLA_Template 
    templ = NULL;

  PLA_Obj_local_width( A, &local_size );
  if ( local_size == 0 ) return;
  PLA_Obj_local_length( A, &local_size );
  if ( local_size == 0 ) return;

  PLA_Obj_view_all( A, &ABR );
  PLA_Obj_view_all( Xdpmv, &Xdpmv_cur );
  PLA_Obj_view_all( Xdpmv_conj_trans, &Xdpmv_conj_trans_cur );

  PLA_Obj_template( A, &templ );
  PLA_Temp_comm_row_rank( templ, &mycol );

  while ( TRUE ) {
    PLA_Obj_split_size ( ABR, PLA_SIDE_TOP,  &size_top,  &dummy );
    PLA_Obj_split_size ( ABR, PLA_SIDE_LEFT, &size_left, &owner_left );

    if ( ( size = min ( size_top, size_left ) ) == 0 ) break;

    PLA_Obj_split_4( ABR, size, size,     &A11, PLA_DUMMY,
		                          &A21, &ABR );

    PLA_Obj_vert_split_2( Xdpmv_conj_trans_cur, size, 
                                 &Xdpmv_conj_trans_left, &Xdpmv_conj_trans_cur );
    PLA_Obj_horz_split_2( Xdpmv_cur, size,   &Xdpmv_top,
                                             &Xdpmv_cur );

    if ( owner_left == mycol ) {
      PLA_Local_herk ( PLA_LOWER_TRIANGULAR, PLA_NO_TRANS,
                       alpha, Xdpmv_top, beta, A11 );

      PLA_Local_gemm ( PLA_NO_TRANS, PLA_NO_TRANS,
                       alpha, Xdpmv_cur, Xdpmv_conj_trans_left, beta, A21 );
    }
  }

  PLA_Obj_free( &ABR );
  PLA_Obj_free( &A11 );
  PLA_Obj_free( &A21 );
  PLA_Obj_free( &Xdpmv_cur );
  PLA_Obj_free( &Xdpmv_top );
  PLA_Obj_free( &Xdpmv_conj_trans_cur );
  PLA_Obj_free( &Xdpmv_conj_trans_left );

  return( PLA_SUCCESS );
}






