/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

static PLA_Obj A_cpy = NULL;

int PLA_Triangular_At_A_enter( int uplo, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS,
    length_A, width_A, 
    objtype;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Triangular_At_A";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Triangular_At_A_enter" );

  old_size_malloced = PLA_Total_size_malloced( );

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if valid uplo */

    if ( !PLA_Valid_uplo_parameter( uplo ) ) {
      PLA_Warning( "illegal uplo" );
      value--;
    }

    /* Check if A is valid matrix */

    if ( A == NULL || !PLA_Valid_object( A ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    PLA_Obj_objtype( A, &objtype );
    if ( objtype != PLA_MATRIX ){
      PLA_Warning( "Invalid objtype for A" );
      value--;
    }      

    /* Check if datatype is currently supported */
    PLA_Obj_datatype( A, &datatype );
    if ( datatype != MPI_FLOAT && datatype != MPI_DOUBLE &&
	 datatype != MPI_COMPLEX && datatype != MPI_DOUBLE_COMPLEX ){
      PLA_Warning( "datatype not yet supported" );
      value--;
    }

    PLA_Obj_global_length( A, &length_A );
    PLA_Obj_global_width( A, &width_A );

    if ( length_A != width_A ){
      PLA_Warning( "Matrix A must be square" );
      value--;
    }      
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Matrix_create_conf_to( A, &A_cpy );

    PLA_Local_copy( A, A_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Triangular_At_A_exit( int method, PLA_Obj A )

{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];


  PLA_Routine_stack_push( "PLA_Triangular_At_A_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj 
      residual = NULL, norm_A = NULL, norm_residual = NULL, A_tmp = NULL,
      one = NULL, minus_one = NULL;
    double
      tol;
    MPI_Datatype
      datatype;
    PLA_Template
      templ;
    
    PLA_Obj_datatype( A, &datatype );
    PLA_Obj_template( A, &templ );

    if ( datatype == MPI_DOUBLE || datatype == MPI_DOUBLE_COMPLEX )
      tol = 0.0000001;
    else
      tol = 0.0001;

    PLA_Matrix_create_conf_to( A, &residual );

    PLA_Local_copy( A, residual );
    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, 
			     residual );
    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, A_cpy );

    PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS, 1, 1, templ,
			 &norm_A );
    PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS, 1, 1, templ,
			 &norm_residual );

    PLA_Matrix_one_norm( A_cpy, norm_A );

    PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

    PLA_Gemm( PLA_NO_TRANS, PLA_TRANS, minus_one, A_cpy, A_cpy, one, 
	       residual );

    PLA_Matrix_one_norm( residual, norm_residual );

    if ( datatype == MPI_DOUBLE ){
      double *norm_A_p, *norm_residual_p;

      PLA_Obj_local_buffer( norm_A, ( void **) &norm_A_p );
      PLA_Obj_local_buffer( norm_residual, ( void ** ) &norm_residual_p );
    
      if ( *norm_residual_p > tol * *norm_A_p ){
	PLA_Warning( "PLA_Triangular_At_A: large relative error encountered" );
	printf( "norm_1(residual) = %le, norm_1(A) = %le\n", 
		*norm_residual_p, *norm_A_p );
	value--;
      }
    }  
    else
      PLA_Warning("residual check not supported for datatype");

    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &A_tmp );
    PLA_Obj_free( &norm_A );
    PLA_Obj_free( &norm_residual );
    PLA_Obj_free( &minus_one );
    PLA_Obj_free( &one );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Triangular_At_A: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


