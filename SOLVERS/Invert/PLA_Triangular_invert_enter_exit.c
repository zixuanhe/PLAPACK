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

int PLA_Triangular_invert_enter( int uplo, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS,
    length_A, width_A, 
    objtype;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Triangular_invert";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Triangular_invert_enter" );

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
    PLA_Mscalar_create_conf_to( A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_cpy );

    PLA_Copy( A, A_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Triangular_invert_exit( int method, PLA_Obj A )

{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];


  PLA_Routine_stack_push( "PLA_Triangular_invert_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj 
      A_tmp = NULL;
    double
      tol, PLA_Local_abs_diff(), PLA_Local_abs_max(), max_A, diff;
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

    PLA_Mscalar_create_conf_to( A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_tmp );

    PLA_Copy( A, A_tmp );

    PLA_Local_tri_inv( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, A_cpy ); 

    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, A_tmp );
    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, A_cpy );

    max_A = PLA_Local_abs_max( A_tmp );

    diff = PLA_Local_abs_diff( A_tmp, A_cpy );

    if ( diff > tol * max_A ){
      PLA_Warning( "PLA_Triangular_invert: large relative error encountered" );
      printf("diff = %le, max_A = %le\n", diff, max_A );
      value--;
    }      

    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &A_tmp );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Triangular_invert: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


