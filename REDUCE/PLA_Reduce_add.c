/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

PLA_Reduce_add( PLA_Obj source, MPI_Op op, PLA_Obj dest )
{
  int 
    objtype,
    length, width, 
    proj_onto, owner_row, owner_col;

  PLA_Obj 
    temp = NULL, one = NULL;

  PLA_Obj_objtype( dest, &objtype );
  PLA_Create_constants_conf_to( dest, NULL, NULL, &one );

  switch ( objtype ){
  case PLA_MVECTOR:
    PLA_Obj_global_width( dest, &width );
    PLA_Mvector_create_conf_to( dest, width, &temp );
    break;

  case PLA_MATRIX:
    PLA_Obj_project_onto( dest, &proj_onto );
    PLA_Matrix_create_conf_to( dest, &temp );
    PLA_Obj_set_orientation( temp, proj_onto );

    break;

  case PLA_MSCALAR:
    PLA_Obj_owner_row( dest, &owner_row );
    PLA_Obj_owner_col( dest, &owner_col );
    PLA_Mscalar_create_conf_to( dest, owner_row, owner_col, &temp );
    break;

  case PLA_PMVECTOR:
    PLA_Obj_project_onto( dest, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_owner_col( dest, &owner_col );
      PLA_Obj_global_width( dest, &width );
      PLA_Pmvector_create_conf_to( dest, proj_onto, owner_col,
				    width, &temp );
    }
    else {
      PLA_Obj_owner_row( dest, &owner_row );
      PLA_Obj_global_length( dest, &length );
      PLA_Pmvector_create_conf_to( dest, proj_onto, owner_row,
				    length, &temp );
    }
    break;
  }

  PLA_Reduce( source, op, temp );

  if ( op == MPI_SUM )
    PLA_Local_axpy( one, temp, dest );
  else
    PLA_Abort( "operation not yet implemented", __LINE__, __FILE__ );
  
  PLA_Obj_free( &temp );
  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}
