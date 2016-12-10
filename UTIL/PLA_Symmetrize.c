/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Symmetrize( int uplo, PLA_Obj A )
{
  int 
    size;
  PLA_Obj
    A_BR = NULL, A_12 = NULL, A_21 = NULL;

  if ( uplo != PLA_LOWER_TRIANGULAR )
    PLA_Abort( "only uplo == PLA_LOWER_TRIANGULAR implemented",
		__LINE__, __FILE__ );

  PLA_Obj_view_all( A, &A_BR );
  while ( TRUE ){
    PLA_Obj_global_length( A_BR, &size );
    if ( size == 0 ) break;
    
    PLA_Obj_split_4( A_BR, 1, 1, PLA_DUMMY, &A_12,
                                  &A_21,      &A_BR );

    PLA_Obj_set_orientation( A_12, PLA_PROJ_ONTO_ROW );
    
    PLA_Copy( A_21, A_12 );
  }

  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_12 );
  PLA_Obj_free( &A_21 );
  
  return PLA_SUCCESS;
}
