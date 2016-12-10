/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/* PLA_Local_gemv.c */

int PLA_Local_gemv( int trans_A, 
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj x,
		    PLA_Obj beta,  PLA_Obj y )
{
  int 
    value = PLA_SUCCESS,
    local_m, local_n, 
    width_x, width_y,
    lda, stride_x, stride_y, proj_onto;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_x, *buf_y,
    *alphabuf;
  char 
    Trans_A[1];

  /* Perform parameter checking */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_gemv_enter( trans_A, alpha, A, x, beta, y );

  /* Scale y.  Notice: BLAS specification does not require y to 
     be scaled if one of the dimensions of A is zero.  We DO want
     to scale in that case */

  PLA_Local_scal( beta, y );

  /* Extract location of local data */

  PLA_Obj_local_buffer( alpha, &alphabuf);
  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( x, &buf_x);
  PLA_Obj_local_buffer( y, &buf_y);

  /* Get matrix dimensions */

  PLA_Obj_local_length( A, &local_m );
  PLA_Obj_local_width ( A, &local_n );
  PLA_Obj_local_ldim ( A, &lda);

  /* Convert PLAPACK constants to strings recognized by BLAS */

  if( PLA_NO_TRANS == trans_A ) 
    Trans_A[0] = 'N';
  else if( PLA_TRANS == trans_A )
    Trans_A[0] = 'T';
  else if( PLA_CONJUGATE_TRANSPOSE == trans_A ) 
    Trans_A[0] = 'C';
  else if(PLA_CONJ == trans_A ) {
    Trans_A[0] = 'N';
    PLA_Conjugate(A);
  }
  
  /* Extract stride of x */

  /*
  PLA_Obj_project_onto( x, &proj_onto );
  PLA_Obj_objtype( x, &objtype );
  if ( proj_onto == PLA_PROJ_ONTO_COL ) 
    stride_x = 1;
  else
    PLA_Obj_local_ldim( x, &stride_x );
  */

  PLA_Obj_local_width ( x, &width_x );
  if ( width_x > 1 )
    PLA_Obj_local_ldim( x, &stride_x );
  else    
    stride_x = 1;

  /* Extract stride of y */

  /*
  PLA_Obj_project_onto( y, &proj_onto );
  if ( proj_onto == PLA_PROJ_ONTO_COL )
    stride_y = 1;
  else
    PLA_Obj_local_ldim( y, &stride_y );
  */

  PLA_Obj_local_width ( y, &width_y );
  if ( width_y > 1 )
    PLA_Obj_local_ldim( y, &stride_y );
  else    
    stride_y = 1;

  /* Depending on datatype, call appropriate BLAS flavor */

  PLA_Obj_datatype(A, &datatype);
  
  if (0 != local_m && 0 != local_n ){
    if( datatype == MPI_DOUBLE ){
      double d_one = 1.0;

      PLA_dgemv( Trans_A, 
		&local_m, &local_n, 
		(double *) alphabuf, (double *) buf_A, &lda, (double *) buf_x, &stride_x, 
		&d_one, (double *) buf_y, &stride_y); 
    }
    else if( datatype == MPI_FLOAT ){
      float f_one = 1.0;

      PLA_sgemv( Trans_A, 
		&local_m, &local_n, 
		(float *) alphabuf, (float *) buf_A, &lda, (float *) buf_x, &stride_x, 
		&f_one, (float *) buf_y, &stride_y); 
    }
    else if ( datatype == MPI_COMPLEX ){
      PLA_COMPLEX c_one = {1.0, 0.0};

      PLA_cgemv( Trans_A, 
		&local_m, &local_n, 
		(PLA_COMPLEX *) alphabuf, (PLA_COMPLEX *) buf_A, &lda, (PLA_COMPLEX *) buf_x, &stride_x, 
		&c_one, (PLA_COMPLEX *) buf_y, &stride_y); 
    }
    else if( datatype == MPI_DOUBLE_COMPLEX ){
      PLA_DOUBLE_COMPLEX z_one = {1.0, 0.0};

      PLA_zgemv( Trans_A, 
		&local_m, &local_n, 
		(PLA_DOUBLE_COMPLEX *) alphabuf, (PLA_DOUBLE_COMPLEX *) buf_A, &lda, (PLA_DOUBLE_COMPLEX *) buf_x, &stride_x, 
		&z_one, (PLA_DOUBLE_COMPLEX *) buf_y, &stride_y); 
    }
  }
  
  if(PLA_CONJ == trans_A) 
    PLA_Conjugate(A);

  /* Check if all went well */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_gemv_exit( trans_A, alpha, A, x, beta, y );

  return value;
}

/* PLA_Local_symv.c */

int PLA_Local_symv( int uplo,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj x,
		    PLA_Obj beta,  PLA_Obj y )
{
  int 
    value = PLA_SUCCESS,
    local_n, 
    lda, stride_x, stride_y, proj_onto;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_x, *buf_y,
    *alphabuf;
  char 
    Uplo[1];

  /* Perform parameter checking */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_symv_enter( uplo, alpha, A, x, beta, y );

  /* Scale y.  Notice: BLAS specification does not require y to 
     be scaled if one of the dimensions of A is zero.  We DO want
     to scale in that case */

  PLA_Local_scal( beta, y );

  /* Extract location of local data */

  PLA_Obj_local_buffer( alpha, &alphabuf);
  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( x, &buf_x);
  PLA_Obj_local_buffer( y, &buf_y);

  /* Get matrix dimensions */

  PLA_Obj_local_length( A, &local_n );
  PLA_Obj_local_ldim ( A, &lda);

  /* Convert PLAPACK constants to strings recognized by BLAS */

  if( PLA_LOWER_TRIANGULAR == uplo ) 
    Uplo[0] = 'L';
  else 
    Uplo[0] = 'U';

  /* Extract stride of x */

  PLA_Obj_project_onto( x, &proj_onto );
  if ( proj_onto == PLA_PROJ_ONTO_COL )
    stride_x = 1;
  else
    PLA_Obj_local_ldim( x, &stride_x );

  /* Extract stride of y */

  PLA_Obj_project_onto( y, &proj_onto );
  if ( proj_onto == PLA_PROJ_ONTO_COL )
    stride_y = 1;
  else
    PLA_Obj_local_ldim( y, &stride_y );

  /* Depending on datatype, call appropriate BLAS flavor */

  PLA_Obj_datatype(A, &datatype);
  
  if ( 0 != local_n ){
    if( datatype == MPI_DOUBLE ){
      double d_one = 1.0;
      
      PLA_dsymv( Uplo, &local_n, 
		alphabuf, buf_A, &lda, buf_x, &stride_x, 
		&d_one, buf_y, &stride_y); 
    }
    else if( datatype == MPI_FLOAT ){
      float f_one = 1.0;
      
      PLA_ssymv( Uplo, &local_n, 
		alphabuf, buf_A, &lda, buf_x, &stride_x, 
		&f_one, buf_y, &stride_y); 
    }
    else if ( datatype == MPI_COMPLEX ){
      float c_one[2] = {1.0,0.0};
      
      PLA_Abort( "PLA_csymv not yet implemented", __LINE__, __FILE__ );
      /*      PLA_csymv( Uplo, &local_n, 
		alphabuf, buf_A, &lda, buf_x, &stride_x, 
		&c_one, buf_y, &stride_y);  */
    }
    else if( datatype == MPI_DOUBLE_COMPLEX ){
      double z_one[2] = {1.0,0.0};
      
      PLA_Abort( "PLA_zsymv not yet implemented", __LINE__, __FILE__ );
      /*
      PLA_zsymv( Uplo, &local_n, 
		  alphabuf, buf_A, &lda, buf_x, &stride_x, 
		  &z_one, buf_y, &stride_y); 
      */
    }
  }
  
  /* Check if all went well */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_symv_exit( uplo, alpha, A, x, beta, y );

  return value;
}

/* PLA_Local_hemv.c */

int PLA_Local_hemv( int uplo,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj x,
		    PLA_Obj beta,  PLA_Obj y )
{
  int 
    value = PLA_SUCCESS,
    local_n, 
    lda, stride_x, stride_y, proj_onto;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_x, *buf_y,
    *alphabuf;
  char 
    Uplo[1];

  /* Perform parameter checking */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_hemv_enter( uplo, alpha, A, x, beta, y );

  /* Scale y.  Notice: BLAS specification does not require y to 
     be scaled if one of the dimensions of A is zero.  We DO want
     to scale in that case */

  PLA_Local_scal( beta, y );

  /* Extract location of local data */

  PLA_Obj_local_buffer( alpha, &alphabuf);
  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( x, &buf_x);
  PLA_Obj_local_buffer( y, &buf_y);

  /* Get matrix dimensions */

  PLA_Obj_local_length( A, &local_n );
  PLA_Obj_local_ldim ( A, &lda);

  /* Convert PLAPACK constants to strings recognized by BLAS */

  if( PLA_LOWER_TRIANGULAR == uplo ) 
    Uplo[0] = 'L';
  else 
    Uplo[0] = 'U';

  /* Extract stride of x */

  PLA_Obj_project_onto( x, &proj_onto );
  if ( proj_onto == PLA_PROJ_ONTO_COL )
    stride_x = 1;
  else
    PLA_Obj_local_ldim( x, &stride_x );

  /* Extract stride of y */

  PLA_Obj_project_onto( y, &proj_onto );
  if ( proj_onto == PLA_PROJ_ONTO_COL )
    stride_y = 1;
  else
    PLA_Obj_local_ldim( y, &stride_y );

  /* Depending on datatype, call appropriate BLAS flavor */

  PLA_Obj_datatype(A, &datatype);
  
  if ( 0 != local_n ){
    if( datatype == MPI_DOUBLE ){
      double d_one = 1.0;
      
      PLA_Abort( "PLA_dhemv not yet implemented", __LINE__, __FILE__ );
      /*
      PLA_dhemv( Uplo, &local_n, 
		alphabuf, buf_A, &lda, buf_x, &stride_x, 
		&d_one, buf_y, &stride_y); 
      */
    }
    else if( datatype == MPI_FLOAT ){
      float f_one = 1.0;
      
      PLA_Abort( "PLA_dhemv not yet implemented", __LINE__, __FILE__ );
      /*
      PLA_shemv( Uplo, &local_n, 
		alphabuf, buf_A, &lda, buf_x, &stride_x, 
		&f_one, buf_y, &stride_y); 
      */
    }
    else if ( datatype == MPI_COMPLEX ){
      float c_one[2] = {1.0,0.0};
      
      PLA_chemv( Uplo, &local_n, 
		alphabuf, buf_A, &lda, buf_x, &stride_x, 
		&c_one, buf_y, &stride_y); 
    }
    else if( datatype == MPI_DOUBLE_COMPLEX ){
      double z_one[2] = {1.0,0.0};
      
      PLA_zhemv( Uplo, &local_n, 
		  alphabuf, buf_A, &lda, buf_x, &stride_x, 
		  &z_one, buf_y, &stride_y); 
    }
  }
  
  /* Check if all went well */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_hemv_exit( uplo, alpha, A, x, beta, y );

  return value;
}

/* PLA_Local_trmv.c */

int  PLA_Local_trmv(int uplo, int trans, int diag, PLA_Obj A, PLA_Obj x)
{
  int 
    value = PLA_SUCCESS,
    local_n,
    lda, stride_x, proj_onto;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_x;

  char 
    Uplo[1], Trans_A[1], Diag[1];

  /* Perform parameter checking */
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_trmv_enter( uplo, trans, diag, A, x );

  /* Extract location of local data */

  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( x, &buf_x);
  
  /* Get matrix dimensions */

  PLA_Obj_local_length( A, &local_n );
  PLA_Obj_local_ldim ( A, &lda);

  /* Convert PLAPACK constants to strings recognized by BLAS */
  
  if ( PLA_LOWER_TRIANGULAR == uplo )
    Uplo[0] = 'L';
  else 
    Uplo[0] = 'U';

  if( PLA_NO_TRANS == trans) 
    Trans_A[0] = 'N';
  else if( PLA_TRANS == trans )
    Trans_A[0] = 'T';
  else if(PLA_CONJ == trans ) {
    Trans_A[0] = 'N';
    PLA_Conjugate(A);
  }
  else /* if( PLA_CONJ_TRANS == trans ) */
    Trans_A[0] = 'C';
  
  if ( PLA_NONUNIT_DIAG == diag )
    Diag[0] = 'N';
  else 
    Diag[0] = 'U';

  /* Extract stride of x */

  PLA_Obj_project_onto( x, &proj_onto );
  if ( proj_onto == PLA_PROJ_ONTO_COL )
    stride_x = 1;
  else
    PLA_Obj_local_ldim( x, &stride_x );

  /* Depending on datatype, call appropriate BLAS flavor */

  PLA_Obj_datatype(A, &datatype);
  
  if (0 != local_n ){
    if( datatype == MPI_DOUBLE )
      PLA_dtrmv( Uplo, Trans_A, Diag,
		 &local_n, buf_A, &lda, buf_x, &stride_x ); 
    else if( datatype == MPI_FLOAT )
      PLA_strmv( Uplo, Trans_A, Diag,
		 &local_n, buf_A, &lda, buf_x, &stride_x ); 
    else if ( datatype == MPI_COMPLEX )
      PLA_ctrmv( Uplo, Trans_A, Diag,
		 &local_n, buf_A, &lda, buf_x, &stride_x ); 
    else if( datatype == MPI_DOUBLE_COMPLEX )
      PLA_ztrmv( Uplo, Trans_A, Diag,
		 &local_n, buf_A, &lda, buf_x, &stride_x ); 
  }
  
  if(PLA_CONJ == trans ) 
    PLA_Conjugate(A);
  
  /* Check if all went well */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_trmv_exit( uplo, trans, diag, A, x );

  return value;
}

/* PLA_Local_trsv.c */

int  PLA_Local_trsv(int uplo, int trans, int diag, PLA_Obj A, PLA_Obj x)
{
  int 
    value = PLA_SUCCESS,
    local_n,
    lda, stride_x, proj_onto;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_x;

  char 
    Uplo[1], Trans_A[1], Diag[1];

  /* Perform parameter checking */
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_trsv_enter( uplo, trans, diag, A, x );

  /* Extract location of local data */

  PLA_Obj_local_buffer( A, &buf_A);
  PLA_Obj_local_buffer( x, &buf_x);
  
  /* Get matrix dimensions */

  PLA_Obj_local_length( A, &local_n );
  PLA_Obj_local_ldim ( A, &lda);

  /* Convert PLAPACK constants to strings recognized by BLAS */
  
  if ( PLA_LOWER_TRIANGULAR == uplo )
    Uplo[0] = 'L';
  else 
    Uplo[0] = 'U';

  if( PLA_NO_TRANS == trans) 
    Trans_A[0] = 'N';
  else if( PLA_TRANS == trans )
    Trans_A[0] = 'T';
  else if(PLA_CONJUGATE_TRANSPOSE == trans ) 
    Trans_A[0] = 'C';
  else if( PLA_CONJUGATE == trans ) {
    Trans_A[0] = 'N';
    PLA_Conjugate(A);
  }
  
  if ( PLA_NONUNIT_DIAG == diag )
    Diag[0] = 'N';
  else 
    Diag[0] = 'U';

  /* Extract stride of x */

  PLA_Obj_project_onto( x, &proj_onto );
  if ( proj_onto == PLA_PROJ_ONTO_COL )
    stride_x = 1;
  else
    PLA_Obj_local_ldim( x, &stride_x );

  /* Depending on datatype, call appropriate BLAS flavor */

  PLA_Obj_datatype(A, &datatype);
  
#if 1
  if ( Uplo[0] != 'L' && Uplo[0] != 'U' )
    PLA_Abort("Value of Uplo got munged!!", __LINE__, __FILE__);
#endif

  if (0 != local_n ){
    if( datatype == MPI_DOUBLE )
      PLA_dtrsv( Uplo, Trans_A, Diag,
		 &local_n, buf_A, &lda, buf_x, &stride_x ); 
    else if( datatype == MPI_FLOAT )
      PLA_strsv( Uplo, Trans_A, Diag,
		 &local_n, buf_A, &lda, buf_x, &stride_x ); 
    else if ( datatype == MPI_COMPLEX )
      PLA_ctrsv( Uplo, Trans_A, Diag,
		 &local_n, buf_A, &lda, buf_x, &stride_x ); 
    else if( datatype == MPI_DOUBLE_COMPLEX )
      PLA_ztrsv( Uplo, Trans_A, Diag,
		 &local_n, buf_A, &lda, buf_x, &stride_x ); 
  }
  
  if(PLA_CONJ == trans ) 
    PLA_Conjugate(A);
  
  /* Check if all went well */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_trsv_exit( uplo, trans, diag, A, x );

  return value;
}


/* PLA_Local_ger.c */

int PLA_Local_ger( PLA_Obj alpha, PLA_Obj x, PLA_Obj y, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS,
    local_m, local_n, local_width,
    lda, stride_x, stride_y, proj_onto;
  
  MPI_Datatype 
    datatype;
  
  void 
    *buf_A, *buf_x, *buf_y,
    *alphabuf;

  /* Perform parameter checking */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_ger_enter( alpha, x, y, A );

  /* Extract location of local data */

  PLA_Obj_local_buffer( alpha, &alphabuf);
  PLA_Obj_local_buffer( x, &buf_x);
  PLA_Obj_local_buffer( y, &buf_y);
  PLA_Obj_local_buffer( A, &buf_A);

  /* Get matrix dimensions */

  PLA_Obj_local_length( A, &local_m );
  PLA_Obj_local_width ( A, &local_n );
  PLA_Obj_local_ldim ( A, &lda);

  /* Extract stride of x */

  PLA_Obj_local_width ( x, &local_width );
  if ( local_width == 1 )
    stride_x = 1;
  else
    PLA_Obj_local_ldim( x, &stride_x );

  /* Extract stride of y */

  PLA_Obj_local_width ( y, &local_width );
  if ( local_width == 1 )
    stride_y = 1;
  else
    PLA_Obj_local_ldim( y, &stride_y );

  /* Depending on datatype, call appropriate BLAS flavor */

  PLA_Obj_datatype(A, &datatype);
  
  if (0 != local_m && 0 != local_n ){
    if( datatype == MPI_DOUBLE )
      PLA_dger( &local_m, &local_n, 
		 alphabuf, buf_x, &stride_x, buf_y, &stride_y,
		 buf_A, &lda ); 
    else if( datatype == MPI_FLOAT )
      PLA_sger( &local_m, &local_n, 
		alphabuf, buf_x, &stride_x, buf_y, &stride_y,
		buf_A, &lda ); 
    else if ( datatype == MPI_COMPLEX )
      PLA_cger( &local_m, &local_n, 
		 alphabuf, buf_x, &stride_x, buf_y, &stride_y, 
		 buf_A, &lda );  
    else if( datatype == MPI_DOUBLE_COMPLEX )
      PLA_zger( &local_m, &local_n, 
		 alphabuf, buf_x, &stride_x, buf_y, &stride_y,
		 buf_A, &lda ); 
  }
  
  /* Check if all went well */

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_ger_exit( alpha, x, y, A );

  return value;
}

