/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Set_triang_to_zero( int uplo, int diag, PLA_Obj A )
{
  int 
    size_top, size_left, size, dummy;

  PLA_Obj
    A_BR = NULL, A_11 = NULL, A_off = NULL;

  PLA_Obj_view_all( A, &A_BR );

  while ( TRUE ){
    PLA_Obj_split_size( A_BR, PLA_SIDE_TOP,  &size_top,  &dummy );
    PLA_Obj_split_size( A_BR, PLA_SIDE_LEFT, &size_left, &dummy );
    size = min( size_top, size_left );
    if ( size == 0 ) break;

    if ( uplo == PLA_UPPER_TRIANGULAR )
      PLA_Obj_split_4( A_BR, size, size, &A_11, PLA_DUMMY,
			                  &A_off, &A_BR );
    else /* uplo == PLA_LOWER_TRIANGULAR */
      PLA_Obj_split_4( A_BR, size, size, &A_11,      &A_off,
			                  PLA_DUMMY, &A_BR );

    PLA_Local_set_triang_to_zero( uplo, diag, A_11 );

    PLA_Obj_set_to_zero( A_off );
  }

  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_off );
      
  return PLA_SUCCESS;
}

int PLA_Local_set_triang_to_zero( int uplo, int diag, PLA_Obj A )
{
  int 
    m, n;

  PLA_Obj
    A_BR = NULL, A_11 = NULL, A_off = NULL;

  PLA_Obj_local_length( A, &m );
  PLA_Obj_local_width ( A, &n );
  if ( m != 0 && n != 0 ) {
    PLA_Obj_view_all( A, &A_BR );

    while ( TRUE ){
      PLA_Obj_local_length( A_BR, &n );
      if ( n == 0 ) break;

      if ( uplo == PLA_UPPER_TRIANGULAR )
	PLA_Obj_split_4( A_BR, 1, 1, &A_11, PLA_DUMMY,
			              &A_off, &A_BR );
      else /* uplo == PLA_LOWER_TRIANGULAR */
	PLA_Obj_split_4( A_BR, 1, 1, &A_11,      &A_off,
			              PLA_DUMMY, &A_BR );

      if ( diag == PLA_UNIT_DIAG ) 
	PLA_Obj_set_to_one( A_11 );

      PLA_Obj_set_to_zero( A_off );
    }
  }

  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_off );
      
  return PLA_SUCCESS;
}
