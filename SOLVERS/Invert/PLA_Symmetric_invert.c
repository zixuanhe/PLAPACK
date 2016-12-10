/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Symmetric_invert( int uplo, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS;


  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Symmetric_invert_enter( uplo, A );

  if ( uplo == PLA_LOWER_TRIANGULAR )
    value = PLA_Symmetric_invert_lower( A );
  else
    PLA_Abort( "Method not supported", __LINE__, __FILE__ );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Symmetric_invert_exit( uplo, A );

  return value;
}



