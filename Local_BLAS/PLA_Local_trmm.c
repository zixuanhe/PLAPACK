/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_trmm( int side, int uplo, int trans_A, int diag,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj B )
{
  int 
    local_m, local_n,
    lda, ldb;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_B, *buf_alpha;

  char 
    Trans_A[1], Side[1], Uplo[1], Diag[1];
  
  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( B, &buf_B);
  
  PLA_Obj_local_buffer( alpha, &buf_alpha);
  
  PLA_Obj_local_length( B, &local_m );
  PLA_Obj_local_width ( B, &local_n );

  Side[0] = ( PLA_SIDE_LEFT == side ? 'L': 'R' );
  Uplo[0] = ( PLA_LOWER_TRIANGULAR == uplo ? 'L': 'U' );
  Diag[0] = ( PLA_UNIT_DIAG == diag ? 'U': 'N' );

  if( PLA_NO_TRANS == trans_A) {
    Trans_A[0] = 'N';
  }
  else if( PLA_TRANS == trans_A ){
    Trans_A[0] = 'T';
  }
  else if(PLA_CONJ == trans_A) {
    Trans_A[0] = 'N';
    PLA_Conjugate(A);
  }
  else /* if( PLA_CONJ_TRANS == trans_A ) */
    {
      Trans_A[0] = 'C';
    }
  
  PLA_Obj_local_ldim ( A, &lda);
  PLA_Obj_local_ldim ( B, &ldb);
  
  PLA_Obj_datatype(A, &datatype);
  
  if (0 != local_m && 0 != local_n ){
    if( datatype == MPI_DOUBLE ){
      double d_one = 1.0;
      
      PLA_dtrmm( Side, Uplo, Trans_A, Diag,
		&local_m, &local_n, 
		buf_alpha, buf_A, &lda, buf_B, &ldb ); 
    }
    else if( datatype == MPI_FLOAT ){
      float f_one = 1.0;
      
      PLA_dtrmm( Side, Uplo, Trans_A, Diag,
		&local_m, &local_n, 
		buf_alpha, buf_A, &lda, buf_B, &ldb ); 
    }
    else if ( datatype == MPI_COMPLEX ){
      float c_one[2] = {1.0,0.0};
      
      PLA_dtrmm( Side, Uplo, Trans_A, Diag,
		&local_m, &local_n, 
		buf_alpha, buf_A, &lda, buf_B, &ldb ); 
    }
    else if( datatype == MPI_DOUBLE_COMPLEX ){
      double z_one[2] = {1.0,0.0};
      
      PLA_dtrmm( Side, Uplo, Trans_A, Diag,
		&local_m, &local_n, 
		buf_alpha, buf_A, &lda, buf_B, &ldb ); 
    }
  }
  
  if(PLA_CONJ == trans_A) 
    PLA_Conjugate(A);
  
  return PLA_SUCCESS;
}
