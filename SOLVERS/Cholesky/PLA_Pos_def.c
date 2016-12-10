/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/**********************************************************
  ROUTINE:
  int PLA_Pos_def_solve ( PLA_Obj A, PLA_Obj B)
 
  PURPOSE:
  Driver routine to compute the solution, X, of A X = B, 
  overwriting B.  A is a positive definite, n x n matrix,
  and B is either a single or multiple right hand side.

  INPUTS:
  A                       PLA_MATRIX, positive definite.
                          On output, the lower triangular 
			  portion of A is overwritten with
			  L.

  B                       PLA_MVECTOR or PLA_PMVECTOR.  

  OUTPUTS:
  A                       The lower triangular portion of A
                          is overwritten with L.

  B                       B is overwritten with the solution.

  ALGORITHM:
  Computes the cholesky factorization of A, A = L L^T,
  then does two triangular solves : Y = L^{-1} B, and
  X = L^{-T} Y.

  
**********************************************************/

int PLA_Pos_def_solve ( PLA_Obj A, PLA_Obj B)
{
  PLA_Obj
    one = NULL;

  int 
    objtype_B,
    trans,
    value = PLA_SUCCESS;

  MPI_Datatype
    datatype_B;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Pos_def_solve_enter( A, B );

  PLA_Obj_objtype( B, &objtype_B );
  PLA_Obj_datatype( B, &datatype_B );
  if ( datatype_B == MPI_DOUBLE || datatype_B == MPI_FLOAT )
    trans = PLA_TRANSPOSE;
  else
    trans = PLA_CONJUGATE_TRANSPOSE;

  PLA_Create_constants_conf_to( A, NULL, NULL, &one);

  value = PLA_Chol ( PLA_LOWER_TRIANGULAR, A );
  
  if (value != PLA_SUCCESS)
      PLA_Warning( "irregular return from PLA_Chol" );

  if ( objtype_B == PLA_MATRIX ){
    PLA_Trsm ( PLA_SIDE_LEFT, PLA_LOWER_TRIANGULAR, 
		PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG,
		one, A, B );

    PLA_Trsm ( PLA_SIDE_LEFT, PLA_LOWER_TRIANGULAR, 
		trans,          PLA_NONUNIT_DIAG,
		one, A, B );
  }
  else if ( objtype_B == PLA_MVECTOR ){
    PLA_Trsv ( PLA_LOWER_TRIANGULAR, 
		PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG,
		A, B );

    PLA_Trsv ( PLA_LOWER_TRIANGULAR, 
		trans, PLA_NONUNIT_DIAG,
		A, B );
  }

  PLA_Obj_free ( &one);

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Pos_def_solve_exit( A, B );

  return value;
}
