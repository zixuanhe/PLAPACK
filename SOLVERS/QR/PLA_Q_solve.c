/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Q_solve( int side, int trans, PLA_Obj A, PLA_Obj s, PLA_Obj B)
/*
  Purpose: Solve one of the following:

       OPERATION             SIDE               TRANS
       ----------------------------------------------
       Q   X = B               PLA_SIDE_LEFT      PLA_NO_TRANS
       Q^T X = B               PLA_SIDE_LEFT      PLA_TRANS
       X Q   = B               PLA_SIDE_RIGHT     PLA_NO_TRANS
       X Q^T = B               PLA_SIDE_RIGHT     PLA_TRANS
       
    where Q is stored as Householder vectors below the diagonal of A and in s
    and B is given.  The solution X overwrites B.

  Input:  side    --   Indicates whether Q appears to left or right of X
                       (int)
          trans   --   Indicates whether to transpose Q
                       (int)
          A       --   General mxn matrix A   
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)
	  B       --   right-hand-side(s)
                       (MVECTOR or MATRIX of length m)

  Output: B       --   first n rows contain least-square solution

  Return value: PLA_SUCCESS iff QR factorization is completed successfully
*/
{
  int       objtype;

  PLA_Obj_objtype ( B, &objtype );
  if ( objtype == PLA_MVECTOR ) 
    return PLA_Q_solve_mv( side, trans, A, s, B );
  else if ( objtype == PLA_MATRIX )
    return PLA_Q_solve_matrix( side, trans, A, s, B );
  else {
    PLA_Abort( "object type unknown in PLA_Q_solve", __LINE__, __FILE__ );
  }
}
