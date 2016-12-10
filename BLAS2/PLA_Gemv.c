/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Gemv( int transa, 
	       PLA_Obj alpha, PLA_Obj A, PLA_Obj x, 
               PLA_Obj beta,  PLA_Obj y )
/*
  Purpose : Parallel matrix vector multiplication

  IN     transa      integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     alpha       multiscalar, scale factor for Ax
  IN     A           matrix
  IN     x           multivector 
  IN     beta        multiscalar, scale factor for y
  IN/OUT y           multivector, overwritten with 
                     y <- alpha * A  * x + beta * y  or 
                     y <- alpha * A' * x + beta * y 

  NOTE:  For details on how to implement parallel matrix-vector multiplication,
         see
	 
	 R. van de Geijn, Using PLAPACK, The MIT Press, 1997.
*/
{
  int 
    value = PLA_SUCCESS,
    owner_row, owner_col;
  PLA_Obj
    alpha_cpy = NULL, beta_cpy = NULL,
    x_dup = NULL, y_dup = NULL, zero = NULL;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Gemv_enter( transa, alpha, A, x, beta, y );

  /* Create usual constants */
  PLA_Create_constants_conf_to( A, NULL, &zero, NULL );

  if ( value == PLA_SUCCESS ){
    /* If necessary, duplicate alpha and beta to all nodes */
    PLA_Obj_owner_row( alpha, &owner_row );
    PLA_Obj_owner_col( alpha, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &alpha_cpy );
      PLA_Copy( alpha, alpha_cpy );
    }

    PLA_Obj_owner_row( beta, &owner_row );
    PLA_Obj_owner_col( beta, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( beta, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &beta_cpy );
      PLA_Copy( beta, beta_cpy );
    }

    if ( transa == PLA_NO_TRANS ){
      /* Create duplicated projected multivectors to hold copy of x and
         local contributions to y */
      PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1,
				    &x_dup );
      PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1,
				    &y_dup );
      
      /* Duplicate x */
      PLA_Copy( x, x_dup );

      /* Perform local matrix-vector multiply */
      PLA_Local_gemv( transa, 
		       ( alpha_cpy == NULL ? alpha: alpha_cpy ), 
		       A, x_dup, zero, y_dup );
	
      /* Add local contributions to vector y */
      PLA_Reduce_x( PLA_SHAPE_GENERAL, y_dup, 
		     ( beta_cpy == NULL ? beta : beta_cpy ), y );
    }
    else{
      /* Create duplicated projected multivectors to hold copy of x and
         local contributions to y */
      PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1,
				    &x_dup );
      PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1,
				    &y_dup );
      
      /* Duplicate x */
      PLA_Copy( x, x_dup );

      /* Perform local matrix-vector multiply */
      PLA_Local_gemv( transa, 
		       ( alpha_cpy == NULL ? alpha: alpha_cpy ), 
		       A, x_dup, zero, y_dup );
	
      /* Add local contributions to vector y */
      PLA_Reduce_x( PLA_SHAPE_GENERAL, y_dup, 
		     ( beta_cpy == NULL ? beta : beta_cpy ), y );
    }
  }

  /* free temporary objects and views */
  PLA_Obj_free( &alpha_cpy );
  PLA_Obj_free( &beta_cpy );
  PLA_Obj_free( &x_dup );
  PLA_Obj_free( &y_dup );
  PLA_Obj_free( &zero );

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Gemv_exit( transa, alpha, A, x, beta, y );

  return value;
}


