/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Inv_scal( PLA_Obj alpha, PLA_Obj x )
{
  int 
    value = PLA_SUCCESS,
    owner_row, owner_col;

  PLA_Obj
    alpha_cpy = NULL;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Inv_scal_enter( alpha, x );

  if ( value == PLA_SUCCESS ){
    /* If necessary, duplicate alpha to all nodes that own part of y */
    PLA_Obj_owner_row( alpha, &owner_row );
    PLA_Obj_owner_col( alpha, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &alpha_cpy );
      PLA_Copy( alpha, alpha_cpy );
    }
      
    PLA_Local_inv_scal( ( alpha_cpy == NULL ? alpha : alpha_cpy ),
			   x );

    PLA_Obj_free( &alpha_cpy );
  }

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Inv_scal_enter( alpha, x );

  return value;
}


