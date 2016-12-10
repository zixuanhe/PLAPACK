/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static PLA_Obj alpha_cpy = NULL, x_cpy = NULL;

static int old_size_malloced;

int PLA_Inv_scal_enter( PLA_Obj alpha, PLA_Obj x )
{
  int
    value = PLA_SUCCESS,
    size, 
    objtype;
  char 
    routine_name[ 35 ] = "PLA_Inv_scal";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Inv_scal_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){
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

    /* Check if x is valid object */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              alpha, PLA_ALL_ROWS, PLA_ALL_COLS, &alpha_cpy );
    PLA_Mscalar_create_conf_to( 
              x, PLA_ALL_ROWS, PLA_ALL_COLS, &x_cpy );

    PLA_Copy( alpha, alpha_cpy );
    PLA_Copy( x, x_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Inv_scal_exit( PLA_Obj alpha, PLA_Obj x )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_x, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();

  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Inv_scal_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj x_temp = NULL;

    PLA_Mscalar_create_conf_to( x, PLA_ALL_ROWS, PLA_ALL_COLS, &x_temp );
    PLA_Copy( x, x_temp );

    max_x = PLA_Local_abs_max( x_cpy );

    PLA_Local_inv_scal( alpha_cpy, x_cpy );

    diff = PLA_Local_abs_diff( x_temp, x_cpy );

    if ( diff > 0.000001 * max_x ){
      PLA_Warning( "PLA_Inv_scal: large relative error encountered" );
      value--;
    }      

    PLA_Obj_free( &alpha_cpy );
    PLA_Obj_free( &x_cpy );

    PLA_Obj_free( &x_temp );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Inv_scal: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Inv_scal_exit" ) != 0 )
    PLA_Warning( "PLA_Inv_scal_exit: history stack corrupted" );

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Inv_scal" ) != 0 )
    PLA_Warning( "PLA_Inv_scal: history stack corrupted" );

  return value;
}

