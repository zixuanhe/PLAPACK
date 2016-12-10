/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

PLA_Copy_sym_tridiag_to_msc( int uplo, PLA_Obj A, PLA_Obj tridiag )
{
  int
    size;
  PLA_Obj
    A_BR = NULL, A_11 = NULL, A_21 = NULL, A_21_1 = NULL,
    temp = NULL,
    diag_1 = NULL, diag_B = NULL,
    subdiag_1 = NULL, subdiag_B = NULL;

  if ( uplo != PLA_LOWER_TRIANGULAR )
    PLA_Abort( "uplo != PLA_LOWER_TRIANGULAR not yet implemented",
		__LINE__, __FILE__ );

  PLA_Mscalar_create_conf_to( tridiag, PLA_ALL_ROWS, PLA_ALL_COLS,
			       &temp );
  PLA_Obj_set_to_zero( temp );
  PLA_Obj_view_all( A, &A_BR );
  PLA_Obj_vert_split_2( temp, 1, &diag_B, &subdiag_B );
  
  while ( TRUE ){
    PLA_Obj_global_length( A_BR, &size );
    if ( size == 0 ) break;

    PLA_Obj_split_4( A_BR, 1, 1, &A_11, PLA_DUMMY,
		                  &A_21, &A_BR );

    PLA_Obj_horz_split_2( diag_B, 1, &diag_1,
                                      &diag_B );

    PLA_Obj_horz_split_2( subdiag_B, 1, &subdiag_1,
                                         &subdiag_B );

    PLA_Local_copy( A_11, diag_1 );
    if ( size > 1 ) {
      PLA_Obj_horz_split_2( A_21, 1,   &A_21_1,
                                        PLA_DUMMY );
      PLA_Local_copy( A_21_1, subdiag_1 );
    }
  }

  PLA_Reduce( temp, MPI_SUM, tridiag );

  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_21 );
  PLA_Obj_free( &A_21_1 );
  PLA_Obj_free( &temp );
  PLA_Obj_free( &diag_1 );
  PLA_Obj_free( &diag_B );
  PLA_Obj_free( &subdiag_1 );
  PLA_Obj_free( &subdiag_B );
  
  return PLA_SUCCESS;
}
