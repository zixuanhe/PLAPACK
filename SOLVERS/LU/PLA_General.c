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
  int PLA_General_solve ( PLA_Obj A, PLA_Obj ipiv, PLA_Obj B)
 
  PURPOSE:
  Driver routine to compute the solution, X, of A X = B, 
  overwriting B.  A is a general square matrix, and B is 
  either a single or multiple right hand side.

  INPUTS:
  A                       PLA_MATRIX, square. 

  ipiv                    PLA_MSCALAR, duplicated on all 
                          nodes.

  B                       PLA_MVECTOR or PLA_PMVECTOR.  On
                          output, B is overwritten with the 
			  solution.

  OUTPUTS:			  
  A                       Overwritten with the factors L and
                          U of the row-permuted matrix, P A.
                          L is unit diagonal and is stored in
			  the strictly lower triangular portion,
			  and U is stored in the upper triangular 
			  portion.

  ipiv                    Contains the pivoting information.  Entry
                          i in ipiv is the integer index of the 
			  pivot row at step i, referenced to 
			  row i. If no pivoting is required, 
			  then ipiv is {0,0,0,...}.

  B                       Overwritten with the solution.

  ALGORITHM:
  Computes the LU factorization with pivoting of A, 

  P A = L U 

  applies the pivots to B, 

  A X = B
  (P A) X = (L U) X = P B

  and performs two triangular solves : 

  U X = L^{-1} (P B) 
  X = U^{-1} (L^{-1} (P B)) .  

  The solution overwrites B the first min(m,n) rows of B.
  
**********************************************************/

int PLA_General_solve ( PLA_Obj A, PLA_Obj ipiv, PLA_Obj B)
{
  PLA_Obj
    one = NULL;

  int
    global_width_B = 0,
    error = 0;

  PLA_Create_constants_conf_to ( A, NULL, NULL, &one);

  PLA_Obj_global_width (B, &global_width_B);

  error = PLA_LU ( A, ipiv );

  if ( error != PLA_SUCCESS )
    {
      printf("PLA_General_solve : LU factorization encountered \
a zero pivot.\nPLA_General_solve exiting abnormally.\n");
      return error;
    }

  if ( global_width_B > 0 )
    {
      PLA_Apply_pivots_to_rows ( B, ipiv );

      if ( global_width_B > 1 )
	{
	  PLA_Trsm (PLA_SIDE_LEFT, PLA_LOWER_TRIANGULAR, 
		    PLA_NO_TRANSPOSE, PLA_UNIT_DIAG,
		    one, A, B);

	  PLA_Trsm (PLA_SIDE_LEFT, PLA_UPPER_TRIANGULAR, 
		    PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG,
		    one, A, B);
	}
      else /* global_width_B == 1 */
	{
	  PLA_Trsv (PLA_LOWER_TRIANGULAR, 
		    PLA_NO_TRANSPOSE, PLA_UNIT_DIAG,
		    A, B);

	  PLA_Trsv (PLA_UPPER_TRIANGULAR, 
		    PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG,
		    A, B);
	}
    }
	  

  PLA_Obj_free ( &one);

  return PLA_SUCCESS;
}
