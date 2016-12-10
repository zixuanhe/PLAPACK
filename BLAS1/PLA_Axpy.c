/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Axpy( PLA_Obj alpha, PLA_Obj x, PLA_Obj y )
{
  int 
    value = PLA_SUCCESS,
    objtype_x, objtype_y, size, 
    align_row_x, align_col_x, align_row_y, align_col_y,
    proj_onto_x, proj_onto_y,
    owner_row, owner_col;

  PLA_Obj
    alpha_cpy = NULL, x_cpy = NULL; 
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Axpy_enter( alpha, x, y );

  PLA_Obj_objtype( y, &objtype_y );
    
  switch( objtype_y ){
  case PLA_MVECTOR:
    PLA_Obj_global_width( y, &size );

    /* If necessary, duplicate alpha to all nodes that own part of y */
    PLA_Obj_owner_row( alpha, &owner_row );
    PLA_Obj_owner_col( alpha, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &alpha_cpy );
      PLA_Copy( alpha, alpha_cpy );
    }
      
    /* If necessary, redistribute data in x */
    PLA_Obj_objtype( x, &objtype_x );

    switch( objtype_x ){
    case PLA_MVECTOR:
      PLA_Obj_global_align_row( y, &align_row_y );
      PLA_Obj_global_align_row( x, &align_row_x );
	
      if ( align_row_x == align_row_y ){
	PLA_Local_axpy( ( alpha_cpy == NULL ? alpha : alpha_cpy ),
			 x, y );
	break;
      }
    default:
      PLA_Mvector_create_conf_to( y, size, &x_cpy );
      PLA_Copy( x, &x_cpy );
      PLA_Local_axpy( ( alpha_cpy == NULL ? alpha : alpha_cpy ),
		         x_cpy, y );
    }	  
    break;
  case PLA_MATRIX:
    /* If necessary, duplicate alpha to all nodes that own part of y */
    PLA_Obj_owner_row( alpha, &owner_row );
    PLA_Obj_owner_col( alpha, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &alpha_cpy );
      PLA_Copy( alpha, alpha_cpy );
    }
      
    /* If necessary, redistribute data in x */
    PLA_Obj_objtype( x, &objtype_x );

    switch( objtype_x ){
    case PLA_MATRIX:
      PLA_Obj_global_align_row( y, &align_row_y );
      PLA_Obj_global_align_row( x, &align_row_x );
      
      PLA_Obj_global_align_col( y, &align_col_y );
      PLA_Obj_global_align_col( x, &align_col_x );
      
      PLA_Obj_project_onto( y, &proj_onto_y );
      PLA_Obj_project_onto( x, &proj_onto_x );
      
      if ( align_row_x == align_row_y && 
	  align_col_x == align_col_y &&
	  proj_onto_x == proj_onto_y ){
	PLA_Local_axpy( ( alpha_cpy == NULL ? alpha : alpha_cpy ),
			x, y );
	break;
      }
    default:
      PLA_Matrix_create_conf_to( y, &x_cpy );
      PLA_Obj_set_orientation( x_cpy, proj_onto_y );
      PLA_Copy( x, &x_cpy );
      
      PLA_Local_axpy( ( alpha_cpy == NULL ? alpha : alpha_cpy ),
		      x_cpy, y );
    }	  
    break;
  default:
    PLA_Abort( "Case not yet implemented", __LINE__, __FILE__ );
  }

  PLA_Obj_free( &alpha_cpy );
  PLA_Obj_free( &x_cpy );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Axpy_exit( alpha, x, y );

  return value;
}


