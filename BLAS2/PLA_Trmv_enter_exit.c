/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static PLA_Obj A_cpy = NULL, x_cpy = NULL;

static int old_size_malloced;

int PLA_Trmv_enter( int uplo, int trans, int diag, 
		    PLA_Obj A, PLA_Obj x )
{
  int
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    objtype, proj_onto;
  char 
    routine_name[ 35 ] = "PLA_Trmv";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Trmv_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){
    /* Check if uplo parameter is valid */
    if ( !PLA_Valid_uplo_parameter( uplo) ){
      PLA_Warning( "Invalid parameter uplo" );
      value--;
    }

    /* Check if transpose parameter is valid */
    if ( !PLA_Valid_trans_parameter( trans ) ){
      PLA_Warning( "Invalid parameter trans" );
      value--;
    }


    /* Check if diag parameter is valid */
    if ( !PLA_Valid_diag_parameter( diag ) ){
      PLA_Warning( "Invalid parameter diag" );
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

    PLA_Obj_global_length( A, &length_A );
    PLA_Obj_global_width ( A, &width_A );

    /* Check if A is square */

    if ( length_A != width_A ){
      PLA_Warning( "A is not square" );
      value--;
    }      

    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    PLA_Obj_objtype( x, &objtype );
    if ( objtype != PLA_MATRIX && objtype != PLA_MVECTOR 
                                && objtype != PLA_PMVECTOR ){
      PLA_Warning( "Invalid objtype for x" );
      value--;
    }      

    PLA_Obj_project_onto( x, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_global_length( x, &length_x );
      PLA_Obj_global_width( x, &width_x );
    }
    else {
      PLA_Obj_global_width( x, &length_x );
      PLA_Obj_global_length( x, &width_x );
    }

    if ( width_x != 1 ){
      PLA_Warning( "x is not of width 1" );
      value--;
    }      
      
    /* Check if dimensions match */

    if ( length_A != length_x ){
      PLA_Warning( "dimension of A does not match length of s" );
      value--;
    }      
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_cpy );
    PLA_Mscalar_create_conf_to( 
              x, PLA_ALL_ROWS, PLA_ALL_COLS, &x_cpy );

    PLA_Copy( A, A_cpy );
    PLA_Copy( x, x_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Trmv_exit( int uplo, int trans, int diag,
	           PLA_Obj A, PLA_Obj x )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, max_x, max_all, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Trmv_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj x_temp = NULL;
    MPI_Datatype
      datatype;

    PLA_Mscalar_create_conf_to( x, PLA_ALL_ROWS, PLA_ALL_COLS, &x_temp );
    PLA_Copy( x, x_temp );

    max_A = PLA_Local_abs_max( A_cpy );
    max_x = PLA_Local_abs_max( x_cpy );

    max_all = max( max_A, max_x );

    PLA_Local_trmv( uplo, trans, diag, A_cpy, x_cpy );

    diff = PLA_Local_abs_diff( x_temp, x_cpy );

    PLA_Obj_datatype( A, &datatype );

    if ( datatype == MPI_DOUBLE || datatype == MPI_DOUBLE_COMPLEX ){
      if ( diff > 0.000001 * max_all ){
	PLA_Warning( "PLA_Trmv: large relative error encountered" );
	printf("diff = %lf tolerance = %lf\n", diff, max_all*0.000001 );
	value--;
      }      
    }
    else{ 
      if ( diff > 0.0001 * max_all ){
	PLA_Warning( "PLA_Trmv: large relative error encountered" );
	printf("diff = %lf tolerance = %lf\n", diff, max_all*0.0001 );
	value--;
      }      
    }

    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &x_cpy );

    PLA_Obj_free( &x_temp );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Trmv: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}

