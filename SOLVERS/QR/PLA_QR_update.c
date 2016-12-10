/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_QR_update( PLA_Obj A, PLA_Obj s )

/*
  Purpose: Compute Householder transform based QR factorization of A

  Input:  A       --   General mxn matrix A   
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)

  Output: A       --   QR factorization.  R is stored in upper-triangular
                       portion of A.  Q is stored in vector form below
                       the diagonal of A, with scaling factors in vector s
          s       --   vector of scalars for Householder transforms

  Return value: PLA_SUCCESS iff QR factorization is completed successfully
*/
{
  int value;

  value = PLA_QR_enter( A, s ); 

  if ( value == PLA_SUCCESS )  
    value = PLA_QR_right_update( A, s );

  if ( value == PLA_SUCCESS ) 
    value = PLA_QR_exit( A, s ); 

  return value;
}


