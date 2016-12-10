/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_symm ( int side, int uplo, 
		      PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
		      PLA_Obj beta,  PLA_Obj C )
{
  int 
    local_m, local_n, 
    lda, ldb, ldc;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_B, *buf_C,
    *alphabuf;
  char 
    Side[1], Uplo[1];
  
  PLA_Local_scal( beta, C );
  
  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( B, &buf_B);
  PLA_Obj_local_buffer( C, &buf_C);
  
  PLA_Obj_local_buffer( alpha, &alphabuf);
  
  PLA_Obj_local_length( C, &local_m );
  PLA_Obj_local_width ( C, &local_n );
  
  if( PLA_SIDE_LEFT == side) 
    Side[0] = 'L';
  else
    Side[0] = 'R';

  if ( PLA_LOWER_TRIANGULAR == uplo )
    Uplo[0] = 'L';
  else
    Uplo[0] = 'U';
  
  PLA_Obj_local_ldim ( A, &lda);
  PLA_Obj_local_ldim ( B, &ldb);
  PLA_Obj_local_ldim ( C, &ldc);
  
  PLA_Obj_datatype(C, &datatype);
  
  if (0 != local_m && 0 != local_n ){
    if( datatype == MPI_DOUBLE ){
      double d_one = 1.0;
      
      PLA_dsymm( Side, Uplo,
		&local_m, &local_n, 
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&d_one, buf_C, &ldc); 
    }
    else if( datatype == MPI_FLOAT ){
      float f_one = 1.0;
      
      PLA_ssymm( Side, Uplo,
		&local_m, &local_n, 
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&f_one, buf_C, &ldc); 
    }
    else if ( datatype == MPI_COMPLEX ){
      PLA_COMPLEX c_one = {1.0,0.0};
      
      PLA_csymm( Side, Uplo,
		&local_m, &local_n, 
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&c_one, buf_C, &ldc); 
    }
    else if( datatype == MPI_DOUBLE_COMPLEX ){
      PLA_DOUBLE_COMPLEX z_one = {1.0,0.0};
      
      PLA_zsymm( Side, Uplo,
		&local_m, &local_n, 
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&z_one, buf_C, &ldc); 
    }
  }
  
  return PLA_SUCCESS;
}
