/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Asum( PLA_Obj x, PLA_Obj alpha )
{
  int 
    value = PLA_SUCCESS;

  PLA_Obj
    alpha_local = NULL;

  if ( PLA_ERROR_CHECKING )    
    value = PLA_Asum_enter( x, alpha );

  PLA_Mscalar_create_conf_to 
    ( alpha, PLA_ALL_ROWS, PLA_ALL_COLS, &alpha_local );

  PLA_Local_asum( x, alpha_local );

  PLA_Reduce( alpha_local, MPI_SUM, alpha );

  PLA_Obj_free( &alpha_local );

  if ( PLA_ERROR_CHECKING )   
    value = PLA_Asum_exit( x, alpha );

  return value;
}


