/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_syr2 ( int uplo, 
		PLA_Obj alpha, PLA_Obj x, PLA_Obj y, PLA_Obj A )
{
  int 
    local_m, local_n, 
    lda, stride_x, stride_y, proj_onto;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_x, *buf_y,
    *alphabuf; 
  char 
    Uplo[1];
  
  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( x, &buf_x);
  PLA_Obj_local_buffer( y, &buf_y);
  
  PLA_Obj_local_buffer( alpha, &alphabuf);
  
  PLA_Obj_local_length( A, &local_m );
  PLA_Obj_local_width ( A, &local_n );

  if( PLA_LOWER_TRIANGULAR == uplo) 
    Uplo[0] = 'L';
  else 
    Uplo[0] = 'U';
  
  PLA_Obj_local_ldim ( A, &lda);

  PLA_Obj_project_onto( x, &proj_onto );
  if ( proj_onto == PLA_PROJ_ONTO_COL )
    stride_x = 1;
  else
    PLA_Obj_local_ldim( x, &stride_x );

  PLA_Obj_project_onto( y, &proj_onto );
  if ( proj_onto == PLA_PROJ_ONTO_COL )
    stride_y = 1;
  else
    PLA_Obj_local_ldim( y, &stride_y );

  PLA_Obj_datatype(A, &datatype);
  
  if (0 != local_m && 0 != local_n ){
    if( datatype == MPI_DOUBLE ){
      PLA_dsyr2( Uplo, &local_m, 
		  (double *) alphabuf, (double *) buf_x, &stride_x, 
		  (double *) buf_y, &stride_y, 
		  (double *) buf_A, &lda); 
    }
    else if( datatype == MPI_FLOAT ){
      PLA_ssyr2( Uplo, &local_m, 
		  (float *) alphabuf, (float *) buf_x, &stride_x, 
		  (float *) buf_y, &stride_y, 
		  (float *) buf_A, &lda); 
    }
    else
      PLA_Abort( "complex not supported", __LINE__, __FILE__ );
    /*    else if ( datatype == MPI_COMPLEX ){
      PLA_csyr2( Uplo, &local_m, 
		  (PLA_COMPLEX *) alphabuf, 
		  (PLA_COMPLEX *) buf_x, &stride_x, 
		  (PLA_COMPLEX *) buf_y, &stride_y, 
		  (PLA_COMPLEX *) buf_A, &lda); 
    }
    else if( datatype == MPI_DOUBLE_COMPLEX ){
      PLA_zsyr2( Uplo, &local_m, 
		  (PLA_DOUBLE_COMPLEX *) alphabuf, 
		  (PLA_DOUBLE_COMPLEX *) buf_x, &stride_x, 
		  (PLA_DOUBLE_COMPLEX *) buf_y, &stride_y, 
		  (PLA_DOUBLE_COMPLEX *) buf_A, &lda); 
    }
    */
  }
  
  return PLA_SUCCESS;
}
