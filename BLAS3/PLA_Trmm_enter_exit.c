/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static PLA_Obj alpha_cpy = NULL, A_cpy = NULL, B_cpy = NULL;

static int old_size_malloced;

int PLA_Trmm_enter( int side, int uplo, int transa, int transb, 
	             PLA_Obj alpha, PLA_Obj A, PLA_Obj B, 
                     PLA_Obj beta,  PLA_Obj C )
{
  int
    value = PLA_SUCCESS,
    size, length_A, width_A, length_B, width_B,
    objtype;
  char 
    routine_name[ 35 ] = "PLA_Trmm";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Trmm_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){
    /* Check if transpose parameters are valid */

    if ( !PLA_Valid_side_parameter( side ) ){
      PLA_Warning( "Invalid parameter side" );
      value--;
    }

    if ( !PLA_Valid_uplo_parameter( uplo ) ){
      PLA_Warning( "Invalid parameter uplo" );
      value--;
    }

    if ( !PLA_Valid_trans_parameter( transa ) ){
      PLA_Warning( "Invalid parameter transa" );
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

    /* Check if B is valid matrix */

    if ( B == NULL || !PLA_Valid_object( B ) ) {
      PLA_Warning( "Invalid object B" );
      value--;
    }

    PLA_Obj_objtype( B, &objtype );
    if ( objtype != PLA_MATRIX ){
      PLA_Warning( "Invalid objtype for B" );
      value--;
    }      

    /* Check if matrix dimensions match */
    PLA_Obj_global_length( A, &length_A );
    PLA_Obj_global_width(  A, &width_A );
    PLA_Obj_global_length( B, &length_B );
    PLA_Obj_global_width(  B, &width_B );

    if ( length_A != width_A ){
      PLA_Warning( "A is not square" );
      value--;
    }      

    if ( side == PLA_SIDE_LEFT ){
      if ( length_A != length_B ){
	PLA_Warning( "length of A does not match length of B" );
	value--;
      }      
    }
    else { /* side == PLA_SIDE_RIGHT */
      if ( length_A != width_B ){
	PLA_Warning( "length of A does not match width of B" );
	value--;
      }      
    }
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              alpha, PLA_ALL_ROWS, PLA_ALL_COLS, &alpha_cpy );
    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_cpy );
    PLA_Mscalar_create_conf_to( 
              B, PLA_ALL_ROWS, PLA_ALL_COLS, &B_cpy );

    PLA_Copy( alpha, alpha_cpy );
    PLA_Copy( A, A_cpy );
    PLA_Copy( B, B_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Trmm_exit( int side, int uplo, int transa, int diag,
	             PLA_Obj alpha, PLA_Obj A, PLA_Obj B )
{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, max_B, max_all, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  PLA_Obj
    one = NULL;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Trmm_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj B_temp = NULL;

    PLA_Mscalar_create_conf_to( B, PLA_ALL_ROWS, PLA_ALL_COLS, &B_temp );
    PLA_Copy( B, B_temp );

    max_A = PLA_Local_abs_max( A_cpy );
    max_B = PLA_Local_abs_max( B_cpy );

    PLA_Local_scal( alpha_cpy, B_cpy );

    max_all = max( max_A, max_B );

    PLA_Create_constants_conf_to( A, NULL, NULL, &one );
    PLA_Local_trmm( side, uplo, transa, diag, alpha_cpy, A_cpy, B_cpy );

    diff = PLA_Local_abs_diff( B_temp, B_cpy );

    if ( diff > 0.000001 * max_all ){
      PLA_Warning( "PLA_Trmm: large relative error encountered" );
      value--;

      PLA_Global_show( "B_temp ", B_temp, "%4.1lf ", " " );
      PLA_Global_show( "B_cpy ", B_cpy, "%4.1lf ", " " );
    }      

    PLA_Obj_free( &alpha_cpy );
    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &B_cpy );

    PLA_Obj_free( &B_temp );
    PLA_Obj_free( &one );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Trmm: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Trmm_exit" ) != 0 )
    PLA_Warning( "PLA_Trmm_exit: history stack corrupted" );

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Trmm" ) != 0 )
    PLA_Warning( "PLA_Trmm: history stack corrupted" );

  return value;
}


