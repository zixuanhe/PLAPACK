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

int PLA_Chol_enter( int uplo, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS,
    length_A, width_A, 
    objtype;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Chol";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Chol_enter" );

  old_size_malloced = PLA_Total_size_malloced( );

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if uplo parameter is valid */
    if ( !PLA_Valid_uplo_parameter( uplo ) ){
      PLA_Warning( "Invalid parameter uplo" );
      value--;
    }

    if ( uplo != PLA_LOWER_TRIANGULAR ){
      PLA_Abort( "Only lower triangular factorization supported",
		  __LINE__, __FILE__ );
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

    /* Check if matrix dimensions match */
    PLA_Obj_global_length( A, &length_A );
    PLA_Obj_global_width(  A, &width_A );

    if ( length_A != width_A ){
      PLA_Warning( "A is not square" );
      value--;
    }      

    /* Check if datatype is currently supported */
    PLA_Obj_datatype( A, &datatype );
    if ( datatype != MPI_FLOAT && datatype != MPI_DOUBLE 
	 && datatype != MPI_COMPLEX && datatype != MPI_DOUBLE_COMPLEX ){
      PLA_Warning( "datatype not yet supported" );
      value--;
    }
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_cpy );

    PLA_Copy( A, A_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Chol_exit( int uplo, PLA_Obj A )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];


  PLA_Routine_stack_push( "PLA_Chol_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj A_temp = NULL;
    MPI_Datatype
      datatype;
    double tol;

    PLA_Obj_datatype( A, &datatype );
    if ( datatype == MPI_DOUBLE || datatype == MPI_DOUBLE_COMPLEX )
      tol = 0.000001;
    else
      tol = 0.0001;

    PLA_Mscalar_create_conf_to( A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_temp );
    PLA_Copy( A, A_temp );

    max_A = PLA_Local_abs_max( A_cpy );

    PLA_Local_chol( uplo, A_cpy );

    diff = PLA_Local_abs_diff( A_temp, A_cpy );

    if ( diff > tol * max_A ){
      PLA_Warning( "PLA_Chol: large relative error encountered" );
      value--;
    }      

    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &A_temp );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Chol: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


