/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Trsm( int side, int uplo, int trans, int diag,
	       PLA_Obj alpha, PLA_Obj A, PLA_Obj B ) 
{
  int 
    value = PLA_SUCCESS,
    nb_alg;
  PLA_Template
    templ = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Trsm_enter( side, uplo, trans, diag, alpha, A, B );

  if ( side == PLA_SIDE_LEFT ) {
    if ( uplo == PLA_LOWER_TRIANGULAR )
      PLA_Trsm_left_lower( trans, diag, alpha, A, B );
    else
      PLA_Trsm_left_upper( trans, diag, alpha, A, B );
  }
  else { /* side == PLA_SIDE_RIGHT */
    if ( uplo == PLA_LOWER_TRIANGULAR )
      PLA_Trsm_right_lower( trans, diag, alpha, A, B );
    else
      PLA_Trsm_right_upper( trans, diag, alpha, A, B );
  }      

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Trsm_exit( side, uplo, trans, diag, alpha, A, B );

  return value;
}
