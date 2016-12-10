/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Obj_set_diagonal( PLA_Obj A, PLA_Obj x )
/*
   PLA_Obj_set_diag

   Purpose:  Create diagonal matrix A with diagonal values equal
   to the values in vector x.

   Input:   x --   PLA_MVECTOR of width 1
                   specifies the values to be placed on the diagonal of A

   OUTPUT:  A --   PLA_MATRIX
*/
{
  int 
    size;
  PLA_Obj 
    A_11 = NULL,  A_BR = NULL,
    x_1  = NULL,  x_B  = NULL;

  PLA_Obj_view_all( A, &A_BR );
  PLA_Obj_view_all( x, &x_B );

  while ( TRUE ){
    PLA_Obj_global_length( x_B, &size );
    if ( size == 0 ) break;
    
    PLA_Obj_split_4( A_BR, 1, 1,   &A_11,      PLA_DUMMY,
		                    PLA_DUMMY, &A_BR );

    PLA_Obj_horz_split_2( x_B, 1,  &x_1,
			            &x_B );

    PLA_Local_copy( x_1, A_11 );
  }

  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &x_1 );
  PLA_Obj_free( &x_B );

  return PLA_SUCCESS;
}
