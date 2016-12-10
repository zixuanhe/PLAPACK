/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/*----------------------------------------------------------------------*/

int PLA_Reduce_x( int shape, PLA_Obj Obj_from, 
                    PLA_Obj alpha, PLA_Obj Obj_to )

/************************************************************************
 
  Reduce contents of Obj_from to Obj_to
  Expert version: allows shape and adding to multiple of target

*************************************************************************/

{
  int 
    owner_row,
    owner_col;

  PLA_Obj
    alpha_copy = NULL;

  if ( shape == PLA_SHAPE_GENERAL ){
    PLA_Obj_owner_row( alpha, &owner_row );
    PLA_Obj_owner_col( alpha, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &alpha_copy );
      PLA_Copy( alpha, alpha_copy );
    }
    
    PLA_Local_scal( ( alpha_copy == NULL ? alpha : alpha_copy ),
		     Obj_to );

    PLA_Reduce_add( Obj_from, MPI_SUM, Obj_to );
  }
  else
    PLA_Abort( "Reduce_x: shape != PLA_SHAPE_GENERAL not yet implemented",
		__LINE__, __FILE__ );

  PLA_Obj_free( &alpha_copy );

  return PLA_SUCCESS;
}
