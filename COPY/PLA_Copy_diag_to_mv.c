/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Copy_diag_to_mv( PLA_Obj from, PLA_Obj to )
{
  int 
    size,
    align_from_row, align_from_col, align_to;

  PLA_Obj
    from_BR = NULL,   from_11 = NULL,
    to_B    = NULL,   to_1    = NULL;

  PLA_Obj_global_align_row( from, &align_from_row );
  PLA_Obj_global_align_col( from, &align_from_col );

  PLA_Obj_global_align_row( to, &align_to );

  if ( align_from_row != align_from_col )
    PLA_Abort( "PLA_Copy_diag_to_mv: currently row and column alignment of matrix must equal",
	       __LINE__, __FILE__ );

  if ( align_from_row != align_to )
    PLA_Abort( "PLA_Copy_diag_to_mv: currently alignment of matrix must equal alignment of mv",
	       __LINE__, __FILE__ );

  PLA_Obj_view_all( from, &from_BR );
  PLA_Obj_view_all( to,   &to_B );

  while ( TRUE ){
    PLA_Obj_global_length( from_BR, &size );
    if ( 0 == size ) break;

    PLA_Obj_split_4( from_BR, 1, 1, &from_11,  PLA_DUMMY,
                                    PLA_DUMMY, &from_BR );

    PLA_Obj_horz_split_2( to_B, 1,  &to_1,
                                    &to_B );

    PLA_Local_copy( from_11, to_1 );
  }

  PLA_Obj_free( &from_BR );
  PLA_Obj_free( &from_11 );
  PLA_Obj_free( &to_B );
  PLA_Obj_free( &to_1 );

  return PLA_SUCCESS;
}
