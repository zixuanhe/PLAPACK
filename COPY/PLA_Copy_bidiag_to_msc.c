/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

PLA_Copy_bidiag_to_msc( int uplo, PLA_Obj A, PLA_Obj bidiag )
{
  int
    size, length, width;
  PLA_Obj
    A_BR = NULL, A_11 = NULL, A_12 = NULL, A_12_1 = NULL,
    temp = NULL,
    diag_1 = NULL, diag_B = NULL,
    superdiag_1 = NULL, superdiag_B = NULL;

  if ( uplo != PLA_UPPER_TRIANGULAR )
    PLA_Abort( "uplo != PLA_UPPER_TRIANGULAR not yet implemented",
		__LINE__, __FILE__ );

  PLA_Mscalar_create_conf_to( bidiag, PLA_ALL_ROWS, PLA_ALL_COLS,
			       &temp );
  PLA_Obj_set_to_zero( temp );
  PLA_Obj_view_all( A, &A_BR );
  PLA_Obj_vert_split_2( temp, 1, &diag_B, &superdiag_B );
  
  while ( TRUE ){
    PLA_Obj_global_length( A_BR, &length );
    PLA_Obj_global_width ( A_BR, &width );
    if ( 0 == ( size = min( length, width ) ) ) break;

    PLA_Obj_split_4( A_BR, 1, 1, &A_11,     &A_12,
		                 PLA_DUMMY, &A_BR );

    PLA_Obj_horz_split_2( diag_B, 1, &diag_1,
                                     &diag_B );

    PLA_Obj_horz_split_2( superdiag_B, 1, &superdiag_1,
                                          &superdiag_B );

    PLA_Local_copy( A_11, diag_1 );
    if ( size > 1 ) {
      PLA_Obj_vert_split_2( A_12, 1,   &A_12_1,PLA_DUMMY );
      PLA_Local_copy( A_12_1, superdiag_1 );
    }
  }

  PLA_Reduce( temp, MPI_SUM, bidiag );

  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_12 );
  PLA_Obj_free( &A_12_1 );
  PLA_Obj_free( &temp );
  PLA_Obj_free( &diag_1 );
  PLA_Obj_free( &diag_B );
  PLA_Obj_free( &superdiag_1 );
  PLA_Obj_free( &superdiag_B );
  
  return PLA_SUCCESS;
}
