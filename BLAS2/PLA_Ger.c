/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Ger(  PLA_Obj alpha, PLA_Obj x, PLA_Obj y, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS,
    owner_row, owner_col;
  PLA_Obj
    alpha_cpy = NULL, 
    x_dup = NULL, y_dup = NULL;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Ger_enter( alpha, x, y, A );

    /* If necessary, duplicate alpha and beta to all nodes */
    PLA_Obj_owner_row( alpha, &owner_row );
    PLA_Obj_owner_col( alpha, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &alpha_cpy );
      PLA_Copy( alpha, alpha_cpy );
    }
    
    PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1,
				  &y_dup );
    PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1,
				  &x_dup );
      
    PLA_Copy( x, x_dup );
    PLA_Copy( y, y_dup );

    PLA_Local_ger( ( alpha_cpy == NULL ? alpha: alpha_cpy ), 
		    x_dup, y_dup, A );
	
    PLA_Obj_free( &alpha_cpy );
    PLA_Obj_free( &x_dup );
    PLA_Obj_free( &y_dup );
    
  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Ger_exit( alpha, x, y, A );

  return value;
}


