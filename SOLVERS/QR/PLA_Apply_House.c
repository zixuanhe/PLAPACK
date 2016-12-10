/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

PLA_Apply_House_transform( int side, PLA_Obj A, PLA_Obj u, PLA_Obj beta )
/*
  PLA_Apply_House_transform

  Purpose: Apply Householder transformation.
           side == PLA_SIDE_LEFT       A <- ( I + beta u u^T ) A
           side == PLA_SIDE_RIGHT      A <- A ( I + beta u u^T )
*/
{
  PLA_Obj
    v = NULL, minus_one = NULL, zero = NULL, one = NULL;
  int 
    proj_onto;

  PLA_Create_constants_conf_to( A, NULL, &zero, &one );

  if ( side == PLA_SIDE_LEFT ){
    PLA_Obj_get_orientation( A, &proj_onto );
    PLA_Obj_set_orientation( A, PLA_PROJ_ONTO_ROW );
    
    PLA_Mvector_create_conf_to( A, 1, &v );

    PLA_Obj_set_orientation( A, proj_onto );

    /* Compute v = - A^T u */
    PLA_Gemv( PLA_TRANSPOSE, one, A, u, zero, v );

    /* Update A <- A + beta u v^T */
    PLA_Ger( beta, u, v, A );
  }
  else if ( side == PLA_SIDE_RIGHT ){
    PLA_Mvector_create_conf_to( A, 1, &v );

    /* Compute v = A u */
    PLA_Gemv( PLA_NO_TRANSPOSE, one, A, u, zero, v );

    /* Update A <- A + beta v u^T */
    PLA_Ger( beta, v, u, A );
  }
  else
    PLA_Abort( "illegal value for side", __LINE__, __FILE__ );

  PLA_Obj_free( &v );
  PLA_Obj_free( &zero );
  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}
