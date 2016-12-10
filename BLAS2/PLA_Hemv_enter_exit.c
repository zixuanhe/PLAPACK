/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static PLA_Obj alpha_cpy = NULL, A_cpy = NULL, x_cpy = NULL, 
                beta_cpy  = NULL, y_cpy = NULL;

static int old_size_malloced;

int PLA_Hemv_enter( int uplo,
	             PLA_Obj alpha, PLA_Obj A, PLA_Obj x, 
                     PLA_Obj beta,  PLA_Obj y )
{
  int
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    objtype, proj_onto;
  char 
    routine_name[ 35 ] = "PLA_Hemv";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Hemv_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){
    /* Check if uplo parameters are valid */
    if ( !PLA_Valid_uplo_parameter( uplo ) ){
      PLA_Warning( "Invalid parameter uplo" );
      value--;
    }

    /* Check if alpha is valid multiscalar of size 1x1 */

    if ( alpha == NULL || !PLA_Valid_object( alpha ) ) {
      PLA_Warning( "Invalid object alpha" );
      value--;
    }

   PLA_Obj_objtype( alpha, &objtype );
    if ( objtype != PLA_MSCALAR ){
      PLA_Warning( "Invalid objtype for alpha" );
      value--;
    }      

    PLA_Obj_global_length( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid global length for alpha" );
      value--;
    }      

    PLA_Obj_global_width( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid global width for alpha" );
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
      
    /* Check if beta is valid multiscalar of size 1x1 */

    if ( beta == NULL || !PLA_Valid_object( beta ) ) {
      PLA_Warning( "Invalid object beta" );
      value--;
    }

    PLA_Obj_objtype( beta, &objtype );
    if ( objtype != PLA_MSCALAR ){
      PLA_Warning( "Invalid objtype for beta" );
      value--;
    }      

    PLA_Obj_global_length( beta, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid global length for beta" );
      value--;
    }      

    PLA_Obj_global_width( beta, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid global width for beta" );
      value--;
    }      

    /* Check if y is valid vector */

    if ( y == NULL || !PLA_Valid_object( y ) ) {
      PLA_Warning( "Invalid object y" );
      value--;
    }

    PLA_Obj_objtype( y, &objtype );
    if ( objtype != PLA_MATRIX && objtype != PLA_MVECTOR 
                                && objtype != PLA_PMVECTOR ){
      PLA_Warning( "Invalid objtype for y" );
      value--;
    }      

    PLA_Obj_project_onto( y, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_global_length( y, &length_y );
      PLA_Obj_global_width( y, &width_y );
    }
    else {
      PLA_Obj_global_width( y, &length_y );
      PLA_Obj_global_length( y, &width_y );
    }

    if ( width_y != 1 ){
      PLA_Warning( "y is not of width 1" );
      value--;
    }      
      
    /* Check if dimensions match */

    if ( length_A != length_y ){
      PLA_Warning( "length of A does not match length of y" );
      value--;
    }      

    if ( width_A != length_x ){
      PLA_Warning( "width of A does not match length of x" );
      value--;
    }      
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              alpha, PLA_ALL_ROWS, PLA_ALL_COLS, &alpha_cpy );
    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_cpy );
    PLA_Mscalar_create_conf_to( 
              x, PLA_ALL_ROWS, PLA_ALL_COLS, &x_cpy );
    PLA_Mscalar_create_conf_to( 
              beta, PLA_ALL_ROWS, PLA_ALL_COLS, &beta_cpy );
    PLA_Mscalar_create_conf_to( 
              y, PLA_ALL_ROWS, PLA_ALL_COLS, &y_cpy );

    PLA_Copy( alpha, alpha_cpy );
    PLA_Copy( A, A_cpy );
    PLA_Copy( x, x_cpy );
    PLA_Copy( beta, beta_cpy );
    if ( PLA_Local_equal_zero( beta_cpy ) )
      PLA_Obj_set_to_zero( y_cpy );
    else
      PLA_Copy( y, y_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Hemv_exit( int uplo,
	             PLA_Obj alpha, PLA_Obj A, PLA_Obj x, 
                     PLA_Obj beta,  PLA_Obj y )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, max_x, max_y, max_all, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  PLA_Obj
    one = NULL;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Hemv_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj y_temp = NULL;
    double 
      tol;
    MPI_Datatype
      datatype;

    PLA_Obj_datatype( A, &datatype );
    if ( datatype == MPI_DOUBLE || datatype == MPI_DOUBLE_COMPLEX )
      tol = 0.000001;
    else
      tol = 0.0001;

    PLA_Mscalar_create_conf_to( y, PLA_ALL_ROWS, PLA_ALL_COLS, &y_temp );
    PLA_Copy( y, y_temp );

    max_A = PLA_Local_abs_max( A_cpy );
    max_x = PLA_Local_abs_max( x_cpy );

    max_y = PLA_Local_abs_max( y_cpy ); 
    max_all = max( max_A, max( max_x, max_y ) ); 

    PLA_Create_constants_conf_to( A, NULL, NULL, &one );
    PLA_Local_hemv( uplo, alpha_cpy, A_cpy, x_cpy,
		           beta_cpy,  y_cpy );

    diff = PLA_Local_abs_diff( y_temp, y_cpy );

    if ( diff > tol * max_all ){
      PLA_Warning( "PLA_Hemv: large relative error encountered" );
      value--;
    }      

    PLA_Obj_free( &alpha_cpy );
    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &x_cpy );
    PLA_Obj_free( &beta_cpy );
    PLA_Obj_free( &y_cpy );

    PLA_Obj_free( &y_temp );
    PLA_Obj_free( &one );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Hemv: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Hemv_exit" ) != 0 )
    PLA_Warning( "PLA_Hemv_exit: history stack corrupted" );

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Hemv" ) != 0 )
    PLA_Warning( "PLA_Hemv: history stack corrupted" );

  return value;
}

