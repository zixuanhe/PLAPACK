/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Syrk ( int uplo, int trans, PLA_Obj alpha, PLA_Obj A,
                                     PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS,
    nb_alg;
  
  PLA_Template
    templ;

  if ( PLA_ERROR_CHECKING )
    value = PLA_Syrk_enter ( uplo, trans, alpha, A, beta, C ); 

  PLA_Obj_template( C, &templ ); 
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  PLA_Syrk_panpan ( nb_alg, uplo, trans, alpha, A, beta, C );

  if ( PLA_ERROR_CHECKING )
    value = PLA_Syrk_exit ( uplo, trans, alpha, A, beta, C ); 

  return PLA_SUCCESS;
}
