/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_gemm( int trans_A, int trans_B,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
		    PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS,
    local_m, local_n, local_k,
    lda, ldb, ldc;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_B, *buf_C,
    *alphabuf;
  char 
    Trans_A[1], Trans_B[1];
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_gemm_enter( trans_A, trans_B, alpha, A, B, beta, C );

  /* Scale matrix C */

  PLA_Local_scal( beta, C );
  
  /* Extract location of local data */

  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( B, &buf_B);
  PLA_Obj_local_buffer( C, &buf_C);
  
  PLA_Obj_local_buffer( alpha, &alphabuf);
  
  /* Get matrix dimensions */

  PLA_Obj_local_length( C, &local_m );
  PLA_Obj_local_width ( C, &local_n );
  
  /* Convert PLAPACK constants to strings recognized by BLAS */

  if( PLA_NO_TRANS == trans_A) {
    Trans_A[0] = 'N';
    PLA_Obj_local_width( A, &local_k );
  }
  else if( PLA_TRANS == trans_A ){
    Trans_A[0] = 'T';
    PLA_Obj_local_length( A, &local_k );
  }
  else if(PLA_CONJ == trans_A) {
    Trans_A[0] = 'N';
    PLA_Obj_local_width( A, &local_k );
    PLA_Conjugate(A);
  }
  else /* if( PLA_CONJ_TRANS == trans_A ) */
    {
      Trans_A[0] = 'C';
      PLA_Obj_local_length(  A, &local_k );
    }
  
  if (PLA_NO_TRANS == trans_B) 
    Trans_B[0] = 'N';
  else if( PLA_TRANS == trans_B)
    Trans_B[0] = 'T';
  else if (PLA_CONJ == trans_B) {
    Trans_B[0] = 'N';
    PLA_Conjugate(B);
  }
  else /* if( PLA_CONJ_TRANS == trans_A ) */
    Trans_B[0] = 'C';
  
  PLA_Obj_local_ldim ( A, &lda);
  PLA_Obj_local_ldim ( B, &ldb);
  PLA_Obj_local_ldim ( C, &ldc);
  
  PLA_Obj_datatype(C, &datatype);
  
  /* Depending on datatype, call appropriate BLAS flavor */

  if (0 != local_m && 0 != local_n && 0 != local_k ){
    if( datatype == MPI_DOUBLE ){
      double d_one = 1.0;

#ifdef PLA_OMP
      PLA_OMP_dgemm(
	        Trans_A, Trans_B,
		local_m, local_n, local_k,
		alphabuf, buf_A, lda, buf_B, ldb, 
	 	&d_one, buf_C, ldc); 
#else      
      PLA_dgemm(
	        Trans_A, Trans_B,
		&local_m, &local_n, &local_k,
		alphabuf, buf_A, &lda, buf_B, &ldb, 
	 	&d_one, buf_C, &ldc); 
#endif
    }
    else if( datatype == MPI_FLOAT ){
      float f_one = 1.0;

#ifdef PLA_OMP
      PLA_OMP_sgemm( Trans_A, Trans_B,
		local_m, local_n, local_k,
		alphabuf, buf_A, lda, buf_B, ldb, 
		&f_one, buf_C, ldc); 
#else
      PLA_sgemm( Trans_A, Trans_B,
		&local_m, &local_n, &local_k,
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&f_one, buf_C, &ldc); 
#endif
    }
    else if ( datatype == MPI_COMPLEX ){
      PLA_COMPLEX c_one = {1.0,0.0};
      
      PLA_cgemm( Trans_A, Trans_B,
		&local_m, &local_n, &local_k,
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&c_one, buf_C, &ldc); 
    }
    else if( datatype == MPI_DOUBLE_COMPLEX ){
      PLA_DOUBLE_COMPLEX z_one = {1.0,0.0};
      
      PLA_zgemm( Trans_A, Trans_B,
		&local_m, &local_n, &local_k,
		alphabuf, buf_A, &lda, buf_B, &ldb, 
		&z_one, buf_C, &ldc); 
    }
  }
  
  if(PLA_CONJ == trans_A) 
    PLA_Conjugate(A);
  
  if (PLA_CONJ == trans_B) 
    PLA_Conjugate(B);
  
  /* Check if all went well */

  if ( PLA_ERROR_CHECKING )
    value = PLA_Local_gemm_exit( trans_A, trans_B, alpha, A, B, beta, C );

  return value;
}
