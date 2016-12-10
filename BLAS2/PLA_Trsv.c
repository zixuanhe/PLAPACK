/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/****************************************************************************/

int PLA_Trsv( int uplo, int trans, int diag, PLA_Obj A, PLA_Obj x)
{
  int 
    value = PLA_SUCCESS,
    objtype;

  PLA_Obj 
    x_temp = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Trsv_enter( uplo, trans, diag, A, x );

  if ( !value ){
    PLA_Obj_objtype( x, &objtype );
    if ( objtype != PLA_MVECTOR ){
      PLA_Mvector_create_conf_to( A, 1, &x_temp );
      PLA_Copy(x, x_temp);
    }

    if ( PLA_UPPER_TRIANGULAR == uplo ){
      if ( PLA_NO_TRANSPOSE == trans )
         PLA_JVTrsv_un( diag, A, ( x_temp == NULL ? x : x_temp ) );  
      else /* PLA_TRANSPOSE */
	/*         PLA_Trsv_ut( diag, A, x_temp ); */
	PLA_Abort( "Upper transpose not yet supported", __LINE__, __FILE__ );
    }
    else /* PLA_LOWER_TRIANGULAR */ {
      if ( PLA_NO_TRANSPOSE == trans )
         PLA_JVTrsv_ln( diag, A, ( x_temp == NULL ? x : x_temp ) );  
      else if ( PLA_TRANSPOSE == trans )
	PLA_JVTrsv_lt( diag, A, ( x_temp == NULL ? x : x_temp ) );  
      else if ( PLA_CONJUGATE_TRANSPOSE == trans )
	PLA_JVTrsv_lc( diag, A, ( x_temp == NULL ? x : x_temp ) );  
    }

    if ( x_temp != NULL )
      PLA_Copy( x_temp, x );

    PLA_Obj_free( &x_temp);
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Trsv_exit( uplo, trans, diag, A, x );

  return value;
}

