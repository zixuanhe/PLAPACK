/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Symv( int uplo,
	       PLA_Obj alpha, PLA_Obj A, PLA_Obj x, 
               PLA_Obj beta,  PLA_Obj y )
{
  int 
    value = PLA_SUCCESS,
    owner_row, owner_col;
  PLA_Obj
    alpha_cpy = NULL, beta_cpy = NULL,
    x_dup_onto_rows = NULL, x_dup_onto_cols = NULL, 
    y_dup_onto_rows = NULL, y_dup_onto_cols = NULL, 
    zero = NULL;

  if ( PLA_ERROR_CHECKING ) 
    value = PLA_Symv_enter( uplo, alpha, A, x, beta, y ); 

  PLA_Create_constants_conf_to( A, NULL, &zero, NULL );

/*  if ( value == PLA_SUCCESS ){ */
  /* If necessary, duplicate alpha and beta to all nodes */
  PLA_Obj_owner_row( alpha, &owner_row );
  PLA_Obj_owner_col( alpha, &owner_col );
  if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
    PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS,
				&alpha_cpy );
    PLA_Copy( beta, beta_cpy );
  }
  PLA_Obj_owner_row( beta, &owner_row );
  PLA_Obj_owner_col( beta, &owner_col );
  if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
    PLA_Mscalar_create_conf_to( beta, PLA_ALL_ROWS, PLA_ALL_COLS,
				&beta_cpy );
    PLA_Copy( beta, beta_cpy );
  }

  PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1,
			        &x_dup_onto_rows );
  PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1,
			        &x_dup_onto_cols );
  PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1,
			        &y_dup_onto_rows );
  PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1,
			        &y_dup_onto_cols );

  PLA_Obj_set_to_zero( y_dup_onto_rows );
  PLA_Obj_set_to_zero( y_dup_onto_cols );

  PLA_Copy( x, x_dup_onto_rows );
  PLA_Copy( x, x_dup_onto_cols );

  PLA_Local_scal( (beta_cpy == NULL ? beta : beta_cpy ), y );

  PLA_Trmv_perform_local_part( uplo, PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG, 
			        A, x_dup_onto_rows, 
			        y_dup_onto_cols );

  PLA_Local_scal( alpha, y_dup_onto_cols );

  PLA_Reduce_add( y_dup_onto_cols, MPI_SUM, y );

  PLA_Trmv_perform_local_part( uplo, PLA_TRANSPOSE, PLA_ZERO_DIAG, 
			        A, x_dup_onto_cols, 
			        y_dup_onto_rows );

  PLA_Local_scal( alpha, y_dup_onto_rows );

  PLA_Reduce_add( y_dup_onto_rows, MPI_SUM, y );

  PLA_Obj_free( &alpha_cpy );
  PLA_Obj_free( &beta_cpy );
  PLA_Obj_free( &x_dup_onto_rows );
  PLA_Obj_free( &x_dup_onto_cols );
  PLA_Obj_free( &y_dup_onto_rows );
  PLA_Obj_free( &y_dup_onto_cols );
  PLA_Obj_free( &zero );

  if ( PLA_ERROR_CHECKING )    
    value = PLA_Symv_exit( uplo, alpha, A, x, beta, y ); 

  return value;
}


