/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_tri_inv( int uplo, int diag, PLA_Obj A )
{
  int 
    local_length, local_width, local_ldim;

  MPI_Datatype
    datatype;

  void 
    *local_buf;

  PLA_Obj_local_length( A, &local_length );
  PLA_Obj_local_width( A, &local_width );
  
  if ( local_length != 0 && local_width != 0 ) {
    PLA_Obj_datatype( A, &datatype );
    PLA_Obj_local_ldim   ( A, &local_ldim );
    PLA_Obj_local_buffer( A, &local_buf );

    if ( datatype == MPI_DOUBLE ){
      if ( uplo == PLA_LOWER_TRIANGULAR )
	PLA_dtri_inv_lower(  diag, local_length, 
			      (double *) local_buf, local_ldim );
      else
	PLA_Abort( "datatype not yet implemented", __LINE__, __FILE__ );
    }
    else
      PLA_Abort( "datatype not yet implemented", __LINE__, __FILE__ );
  }

  return PLA_SUCCESS;
}

#define A( i,j,lda )  a[ ((j)-1) * (lda) + (i)-1 ]

PLA_dtri_inv_lower( int diag, int n, double *a, int lda )
{
  int 
    n_1, n_2;

  double
    d_minus_one = -1.0, d_one = 1.0;

  if ( n == 1 )
    A( 1, 1, lda ) = 1.0/ A( 1,1,lda );
  else{
    n_1 = n/2;
    n_2 = n - n_1;

    PLA_dtrsm( "Left", "Lower", "No trans", 
		( diag == PLA_NONUNIT_DIAG ? "Nonunit" : "Unit" ),
		&n_2, &n_1, &d_minus_one, &A( n_1+1,n_1+1,lda ), &lda,
		&A( n_1+1,1,lda ), &lda );

    PLA_dtrsm( "Right", "Lower", "No trans", 
		( diag == PLA_NONUNIT_DIAG ? "Nonunit" : "Unit" ),
		&n_2, &n_1, &d_one, &A( 1,1,lda ), &lda,
		&A( n_1+1,1,lda ), &lda );

    PLA_dtri_inv_lower( diag, n_2, &A( n_1+1,n_1+1,lda ), lda );

    PLA_dtri_inv_lower( diag, n_1, &A( 1, 1, lda ), lda );
  }
  
  return PLA_SUCCESS;
}
    
