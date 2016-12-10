/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Least_sqr_solve( PLA_Obj A, PLA_Obj s, PLA_Obj Q, PLA_Obj B)
/*
  Purpose: Solve linear least squares problem 
                       min || A X - B ||_2 
  overwriting the first n rows of B with X.  The method computes 
  a Householder QR factorization of A, storing the Householder
  vectors below the diagonal in A and in s, updates B <- Q^T B, 
  and solves B <- inv(R) B.  If Q is not equal to PLA_DUMMY,
  Q is computed in that object.

  Input:  A       --   General mxn matrix A   
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)
	  Q       --   General mxm matrix for storing Q if desired.
                       (PLA_MATRIX)
	  B       --   Multivector of length m or general matrix with
                       m rows that stores the right-hand-side(s)
                       ( PLA_MVECTOR or PLA_MATRIX )
  Output: A       --   QR factorization.  R is stored in upper-triangular
                       portion of A.  Q is stored in vector form below
                       the diagonal of A, with scaling factors in vector s
          s       --   Vector of scalars for Householder transforms
	  Q       --   Matrix Q if Q != PLA_DUMMY
	  B       --   Solution X

  Return value: PLA_SUCCESS iff QR factorization is completed successfully
*/
{
  int     	objtype, width;
  PLA_Obj 	C = NULL, 
                R = NULL, 
		one = NULL;  

  /* Create the usual constants used by the BLAS */
  PLA_Create_constants_conf_to( A, NULL, NULL, &one);

  /* Compute the QR factorization of A, overwriting A with R
     and the Householder vectors.  Scaling factors beta for
     the Householder transforms are stored in s */

  PLA_QR( A, s);

  PLA_Obj_global_width( A, &width );

  /* Determine if B is a single right-hand-side stored in a MVECTOR
     or a matrix */
  PLA_Obj_objtype( B, &objtype );

  /* solve Q C = B overwriting the first width rows of B with C */

  PLA_Q_solve( PLA_SIDE_LEFT, PLA_NO_TRANS, A, s, B);

  /* solve R X = C overwriting C with X */

  PLA_Obj_horz_split_2(A, width, &R,
		                 PLA_DUMMY);
  PLA_Obj_horz_split_2(B, width, &C,
	                         PLA_DUMMY);

  if ( objtype == PLA_MVECTOR )
    PLA_Trsv( PLA_UPPER_TRIANGULAR, PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG, 
	      R, C );
  else
    PLA_Trsm(PLA_SIDE_LEFT,PLA_UPPER_TRIANGULAR, 
             PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG, 
             one, R, C );

  /* If desired, form Q from the Householder vectors stored in A and s */

  if ( Q != PLA_DUMMY ) PLA_Form_Q( PLA_NO_TRANSPOSE, A, s, Q ); 

  PLA_Obj_free( &C );
  PLA_Obj_free( &R );
  PLA_Obj_free(&one);

  return PLA_SUCCESS;
}
