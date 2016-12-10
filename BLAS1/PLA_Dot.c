/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Dot( PLA_Obj x, PLA_Obj y, PLA_Obj alpha )
{
  int 
    value = PLA_SUCCESS;

  PLA_Obj
    alpha_local = NULL;

  if ( PLA_ERROR_CHECKING )    
    value = PLA_Dot_enter( x, y, alpha );

  PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS, 
		               &alpha_local );

  PLA_Local_dot( x, y, alpha_local );

  PLA_Reduce( alpha_local, MPI_SUM, alpha );

  PLA_Obj_free( &alpha_local );

  if ( PLA_ERROR_CHECKING )   
    value = PLA_Dot_exit( x, y, alpha );

  return value;
}


