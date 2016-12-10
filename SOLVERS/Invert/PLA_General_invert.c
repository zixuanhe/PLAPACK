/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_General_invert( int method, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_General_invert_enter( method, A );

  if ( method == PLA_METHOD_INV )
    value = PLA_General_invert_inv( A );
  else if ( method == PLA_METHOD_FACTORS )
    value = PLA_General_invert_factors( A ); 
  else if ( method == PLA_METHOD_STABLE )
    value = PLA_General_invert_stable( A ); 
  else if ( method == PLA_METHOD_TRSM ){
    PLA_Obj
      Temp = NULL, pivots = NULL, one = NULL;
    PLA_Template
      templ;
    int
      size,
      align_row;

    PLA_Create_constants_conf_to( A, NULL, NULL, &one );

    PLA_Matrix_create_conf_to( A, &Temp );
    PLA_Obj_set_to_identity( Temp );

    PLA_Obj_template( A, &templ );
    PLA_Obj_global_length( A, &size );
    PLA_Obj_global_align_row( A, &align_row );
    PLA_Mvector_create( MPI_INT, size, 1, templ, align_row, &pivots );

    PLA_LU( A, pivots );

    PLA_Apply_pivots_to_rows ( Temp, pivots);

    PLA_Trsm( PLA_SIDE_LEFT, PLA_LOWER_TRIANGULAR,
	       PLA_NO_TRANSPOSE, PLA_UNIT_DIAG, one, A, Temp );

    PLA_Trsm( PLA_SIDE_LEFT, PLA_UPPER_TRIANGULAR,
	       PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG, one, A, Temp );
    
    PLA_Local_copy( Temp, A );

    PLA_Obj_free( &Temp );
    PLA_Obj_free( &pivots );
    PLA_Obj_free( &one );

    value = PLA_SUCCESS;
  }
  else
    PLA_Abort( "Method not supported", __LINE__, __FILE__ );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_General_invert_exit( method, A );

  return value;
}



