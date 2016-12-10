/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_LU_right_mv( PLA_Obj A_mv, PLA_Obj A_11_msc, PLA_Obj pivots_msc )
{
  int       
    value = PLA_SUCCESS,
    size, length, width;

  PLA_Obj  
    a_21         = NULL,    A_BR = NULL,     A_B = NULL,
    alpha_11_dup = NULL,    a_12_dup = NULL, a_1 = NULL,
    A_11_BR      = NULL,    A_11_B = NULL,   a_11_1 = NULL,
    pivots_B        = NULL,     pivots_1       = NULL,
    one             = NULL,     minus_one      = NULL;

  PLA_Template
    templ;

  MPI_Datatype
    datatype;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_LU_mv_enter( A_mv, A_11_msc, pivots_msc );

  /* Create usual duplicated scalar constants */
  PLA_Create_constants_conf_to( A_mv, &minus_one, NULL, &one );

  /* View ABR = A */
  PLA_Obj_view_all( A_mv, &A_BR );
  PLA_Obj_view_all( A_mv, &A_B );
  PLA_Obj_view_all( pivots_msc, &pivots_B );

  PLA_Obj_datatype( A_mv, &datatype );
  PLA_Obj_template( A_mv, &templ );
  PLA_Obj_global_length( A_mv, &length );
  PLA_Obj_global_width( A_mv, &width );
  size = min( length, width );

  PLA_Obj_view_all( A_11_msc, &A_11_B );
  PLA_Obj_view_all( A_11_msc, &A_11_BR );

  while ( TRUE ) {
    /* Determine size of current panel */
    PLA_Obj_global_length( A_BR, &length );
    PLA_Obj_global_width( A_BR,  &width );
    if ( min( length, width ) == 0 ) break;

    PLA_Obj_vert_split_2( A_BR, 1, &a_1, PLA_DUMMY );

    PLA_Obj_split_4( A_11_BR, 1, 1,   &alpha_11_dup, &a_12_dup,
                                       PLA_DUMMY,     &A_11_BR );

    PLA_Obj_horz_split_2( pivots_B, 1,    &pivots_1, 
                                           &pivots_B );

    PLA_Iamax( a_1, pivots_1, alpha_11_dup );

    PLA_Obj_horz_split_2( A_11_B, 1, &a_11_1,
                                      &A_11_B );

    Pivot_rows_mv( A_B, pivots_1, a_11_1 );

    PLA_Obj_split_4( A_BR, 1, 1,   PLA_DUMMY, PLA_DUMMY,
                                    &a_21,     &A_BR );

    PLA_Local_inv_scal( alpha_11_dup, a_21 );

    PLA_Local_ger( minus_one, a_21, a_12_dup, A_BR );
    
    PLA_Obj_horz_split_2( A_B, 1,    PLA_DUMMY,
			              &A_B );
  }

  PLA_Obj_free( &a_21 );
  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_B );
  PLA_Obj_free( &alpha_11_dup );
  PLA_Obj_free( &a_12_dup );
  PLA_Obj_free( &a_1 );
  PLA_Obj_free( &A_11_BR );
  PLA_Obj_free( &A_11_B );
  PLA_Obj_free( &a_11_1 );
  PLA_Obj_free( &pivots_B );
  PLA_Obj_free( &pivots_1 );
  PLA_Obj_free( &one );
  PLA_Obj_free( &minus_one );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_LU_mv_exit( A_mv, A_11_msc, pivots_msc );

  return value;
}

