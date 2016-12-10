/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Trmv ( int uplo, int trans, int diag,
		PLA_Obj A, PLA_Obj x )
/*
  Purpose : Parallel triangular matrix vector multiplication

  IN     uplo        integer, PLA_LOWER_TRIANGULAR or PLA_UPPER_TRIANGULAR
  IN     trans       integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     diag        integer, PLA_UNIT_DIAG or PLA_NONUNIT_DIAG
  IN     A           matrix
  IN/OUT x           multivector, overwritten with 
                     x <- A  * x or
                     x <- A' * x

  NOTE:  For details on how to implement parallel matrix-vector multiplication,
         see
	 
	 R. van de Geijn, Using PLAPACK, The MIT Press, 1997.
*/
{
  int 
    value = PLA_SUCCESS;
  PLA_Obj
    x_dup = NULL, y_dup = NULL;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Trmv_enter( uplo, trans, diag, A, x );

  if ( trans == PLA_NO_TRANS ){
    /* Create duplicated projected multivectors to hold copy of x and
       local contributions to y */
    PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1,
				 &x_dup );
    PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1,
				 &y_dup );
      
    /* Duplicate x */
    PLA_Copy( x, x_dup );

    /* Perform local part of triangular matrix-vector multiply */
    PLA_Trmv_perform_local_part( uplo, trans, diag, A, x_dup, y_dup );
	
    /* Add local contributions to vector y */
    PLA_Reduce( y_dup, MPI_SUM, x );
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

    /* Perform local part of triangular matrix-vector multiply */
    PLA_Trmv_perform_local_part( uplo, trans, diag, A, x_dup, y_dup );

    /* Add local contributions to vector x */
    PLA_Reduce( y_dup, MPI_SUM, x );
  }

  /* free temporary objects and views */
  PLA_Obj_free( &x_dup );
  PLA_Obj_free( &y_dup );

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Trmv_exit( uplo, trans, diag, A, x );

  return value;
}


