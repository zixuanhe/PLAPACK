/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Obj_set_to_identity( PLA_Obj A )
{
  int length, width;
  PLA_Obj A11 = NULL,     ABR = NULL;

  PLA_Obj_global_length( A, &length );
  PLA_Obj_global_width( A, &width );

  if ( length != width ){
    printf("PLA_Obj_set_to_identity() error:  object not square\n");
    exit ( 0 );
  }

  PLA_Obj_set_to_zero( A );
  PLA_Obj_view_all( A, &ABR );
  
  while ( TRUE ){
    PLA_Obj_global_length( ABR, &length );
    if ( 0 == length ) break;
    
    PLA_Obj_split_4( ABR, 1, 1, &A11,      PLA_DUMMY,
                                 PLA_DUMMY, &ABR );

    PLA_Obj_set_to_one( A11 );
  }

  PLA_Obj_free( &A11 );
  PLA_Obj_free( &ABR );

  return PLA_SUCCESS;
}
