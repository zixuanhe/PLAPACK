/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

static double Local_BLAS3_start_time;

/***************************************************************************/

int PLA_Local_gemm_enter( int trans_A, int trans_B,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
		    PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS,
    size, 
    length_A, width_A, 
    length_B, width_B, 
    length_C, width_C;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_gemm";

  if ( PLA_TIMINGS ){
    MPI_Barrier( MPI_COMM_WORLD );

    Local_BLAS3_start_time = MPI_Wtime();
  }

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_gemm_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){
    /* Check if trans_A parameter is valid */
    if ( !PLA_Valid_trans_parameter( trans_A ) ){
      PLA_Warning( "Invalid parameter trans_A" );
      value--;
    }

    /* Check if trans_B parameter is valid */
    if ( !PLA_Valid_trans_parameter( trans_B ) ){
      PLA_Warning( "Invalid parameter trans_B" );
      value--;
    }

    /* Check if alpha is valid multiscalar of size 1x1 */

    if ( alpha == NULL || !PLA_Valid_object( alpha ) ) {
      PLA_Warning( "Invalid object alpha" );
      value--;
    } else {

      PLA_Obj_local_length( alpha, &size );
      if ( size != 1 ){
	PLA_Warning( "Invalid local length for alpha" );
	value--;
      }      

      PLA_Obj_local_width( alpha, &size );
      if ( size != 1 ){
	PLA_Warning( "Invalid local width for alpha" );
	value--;
      }
    }

    /* Check if beta is valid multiscalar of size 1x1 */

    if ( beta == NULL || !PLA_Valid_object( beta ) ) {
      PLA_Warning( "Invalid object beta" );
      value--;
    } else {

      PLA_Obj_local_length( beta, &size );
      if ( size != 1 ){
	PLA_Warning( "Invalid local length for beta" );
	value--;
      }      

      PLA_Obj_local_width( beta, &size );
      if ( size != 1 ){
	PLA_Warning( "Invalid local width for beta" );
	value--;
      }      

    }

    /* Check if A is valid matrix */

    if ( A == NULL || !PLA_Valid_object( A ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    /* Check if B is valid matrix */

    if ( B == NULL || !PLA_Valid_object( B ) ) {
      PLA_Warning( "Invalid object B" );
      value--;
    }

    /* Check if C is valid matrix */

    if ( C == NULL || !PLA_Valid_object( C ) ) {
      PLA_Warning( "Invalid object C" );
      value--;
    }

    /* If all parameters are valid, check that they are consistent */

    if ( value == PLA_SUCCESS ) {

      /* Check if matrix dimensions match */

      PLA_Obj_local_length( A, &length_A );
      PLA_Obj_local_width(  A, &width_A );
      PLA_Obj_local_length( B, &length_B );
      PLA_Obj_local_width(  B, &width_B );
      PLA_Obj_local_length( C, &length_C );
      PLA_Obj_local_width(  C, &width_C );

      if ( length_C != (trans_A == PLA_NO_TRANSPOSE ? length_A : width_A) ) {
	PLA_Warning( "First dimension of C does not match first dimension of A" );
	value--;
      }

      if ( width_C != (trans_B == PLA_NO_TRANSPOSE ? width_B : length_B) ) {
	PLA_Warning( "Second dimension of C does not match second dimension of B" );
	value--;
      }

      if (    (trans_A == PLA_NO_TRANSPOSE ?  width_A : length_A) 
	   != (trans_B == PLA_NO_TRANSPOSE ? length_B : width_B) ) {
	PLA_Warning( "Second dimension of A does not match first dimension of B" );
	value--;
      }

    }
 
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_gemm_exit(  int trans_A, int trans_B,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
		    PLA_Obj beta,  PLA_Obj C )
{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_gemm_exit" );

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_GEMM_TIMING] += MPI_Wtime() - Local_BLAS3_start_time;

    Local_BLAS3_start_time = MPI_Wtime();

    MPI_Barrier( MPI_COMM_WORLD );

    PLA_TIMINGS[PLA_LOCAL_GEMM_TIMING_IMBALANCE] += MPI_Wtime() - Local_BLAS3_start_time;
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Local_herk_enter( int uplo, int trans,
		     PLA_Obj alpha, PLA_Obj A, PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    Local_BLAS3_start_time = MPI_Wtime();
  }

  return value;
}


int PLA_Local_herk_exit( int uplo, int trans,
		     PLA_Obj alpha, PLA_Obj A, PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_HERK_TIMING] += MPI_Wtime() - Local_BLAS3_start_time;
  }

  return value;
}

/***************************************************************************/

int PLA_Local_symm_enter( int side, int uplo, 
		      PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
		      PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    Local_BLAS3_start_time = MPI_Wtime();
  }

  return value;
}


int PLA_Local_symm_exit( int side, int uplo, 
		      PLA_Obj alpha, PLA_Obj A, PLA_Obj B,
		      PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_SYMM_TIMING] += MPI_Wtime() - Local_BLAS3_start_time;
  }

  return value;
}

/***************************************************************************/

int PLA_Local_syrk_enter( int uplo, int trans,
		     PLA_Obj alpha, PLA_Obj A, PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    Local_BLAS3_start_time = MPI_Wtime();
  }

  return value;
}

int PLA_Local_syrk_exit( int uplo, int trans,
		     PLA_Obj alpha, PLA_Obj A, PLA_Obj beta,  PLA_Obj C )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_SYRK_TIMING] += MPI_Wtime() - Local_BLAS3_start_time;
  }

  return value;
}

/***************************************************************************/

int PLA_Local_trmm_enter( int side, int uplo, int trans_A, int diag,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj B )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    Local_BLAS3_start_time = MPI_Wtime();
  }

  return value;
}

int PLA_Local_trmm_exit( int side, int uplo, int trans_A, int diag,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj B )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_TRMM_TIMING] += MPI_Wtime() - Local_BLAS3_start_time;
  }

  return value;
}

/***************************************************************************/


int PLA_Local_trsm_enter( int side, int uplo, int trans_A, int diag,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj B )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    Local_BLAS3_start_time = MPI_Wtime();
  }
 
  return value;
}

int PLA_Local_trsm_exit( int side, int uplo, int trans_A, int diag,
		    PLA_Obj alpha, PLA_Obj A, PLA_Obj B )
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_TRSM_TIMING] += MPI_Wtime() - Local_BLAS3_start_time;
  }

  return value;
}

