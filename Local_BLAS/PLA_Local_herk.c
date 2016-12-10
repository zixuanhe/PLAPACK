/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_herk( int uplo, int trans,
		     PLA_Obj alpha, PLA_Obj A, PLA_Obj beta,  PLA_Obj C )
{
  int 
    local_n, local_k,
    lda, ldc;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_C,
    *alphabuf;
  char 
    Trans[1], Uplo[1];
  
  PLA_Local_scal( beta, C );

  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( C, &buf_C);
  
  PLA_Obj_local_buffer( alpha, &alphabuf);
  
  PLA_Obj_local_length( C, &local_n );
  
  if ( PLA_LOWER_TRIANGULAR == uplo )
    Uplo[0] = 'L';
  else if ( PLA_UPPER_TRIANGULAR == uplo )
    Uplo[0] = 'U';

  if( PLA_NO_TRANS == trans) {
    Trans[0] = 'N';
    PLA_Obj_local_width( A, &local_k );
  }
  else if( PLA_TRANS == trans ){
    Trans[0] = 'T';
    PLA_Obj_local_length( A, &local_k );
  }
  else if(PLA_CONJ == trans) {
    Trans[0] = 'N';
    PLA_Obj_local_width( A, &local_k );
    PLA_Conjugate(A);
  }
  else /* if( PLA_CONJ_TRANS == trans ) */
    {
      Trans[0] = 'C';
      PLA_Obj_local_length(  A, &local_k );
    }
  
  PLA_Obj_local_ldim ( A, &lda);
  PLA_Obj_local_ldim ( C, &ldc);
  
  PLA_Obj_datatype(C, &datatype);
  
  if ( 0 != local_n && 0 != local_k ){
    if( datatype == MPI_DOUBLE ){
      double d_one = 1.0;

      PLA_Abort( "PLA_dherk not yet implemented", __LINE__, __FILE__ );
      /*
      PLA_dherk( Uplo, Trans,
		&local_n, &local_k,
		alphabuf, buf_A, &lda, 
		&d_one, buf_C, &ldc); 
      */
    }
    else if( datatype == MPI_FLOAT ){
      float f_one = 1.0;
      
      PLA_Abort( "PLA_dherk not yet implemented", __LINE__, __FILE__ );
      /*
      PLA_sherk( Uplo, Trans,
		&local_n, &local_k,
		alphabuf, buf_A, &lda, 
		&f_one, buf_C, &ldc); 
      */
    }
    else if ( datatype == MPI_COMPLEX ){
      PLA_COMPLEX c_one = {1.0,0.0};
      
      PLA_cherk( Uplo, Trans,
		&local_n, &local_k,
		alphabuf, buf_A, &lda, 
		&c_one, buf_C, &ldc); 
    }
    else if( datatype == MPI_DOUBLE_COMPLEX ){
      PLA_DOUBLE_COMPLEX z_one = {1.0,0.0};
      
      PLA_zherk( Uplo, Trans,
		&local_n, &local_k,
		alphabuf, buf_A, &lda, 
		&z_one, buf_C, &ldc); 
    }
  }
  
  if(PLA_CONJ == trans) 
    PLA_Conjugate(A);
  
  return PLA_SUCCESS;
}
