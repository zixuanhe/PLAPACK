/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include <stdio.h> 

static PLA_Obj x_cpy = NULL, k_cpy = NULL, xmax_cpy = NULL;

static int old_size_malloced;

static double Iamax_start_time;

int PLA_Iamax_enter( PLA_Obj x, PLA_Obj k, PLA_Obj xmax )
{
  int
    value = PLA_SUCCESS,
    size, length_x, width_x, length_k, width_k, length_xmax, width_xmax,
    objtype;
  char 
    routine_name[ 35 ] = "PLA_Iamax";
  MPI_Datatype
    datatype;

  if ( PLA_TIMINGS ){
    Iamax_start_time = MPI_Wtime();
  }

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Iamax_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){
    /* Check if k is valid multiscalar of size 1x1 */

    if ( k == NULL || !PLA_Valid_object( k ) ) {
      PLA_Warning( "Invalid object k" );
      value--;
    }

   PLA_Obj_objtype( k, &objtype );
    if ( objtype != PLA_MSCALAR ){
      PLA_Warning( "Invalid objtype for k" );
      value--;
    }      

   PLA_Obj_datatype( k, &datatype );
    if ( datatype != MPI_INT ){
      PLA_Warning( "Invalid datatype for k" );
      value--;
    }      

    PLA_Obj_global_length( k, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid global length for k" );
      value--;
    }      

    PLA_Obj_global_width( k, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid global width for k" );
      value--;
    }      

    /* Check if xmax is valid multiscalar of size 1x1 */

    if ( xmax == NULL || !PLA_Valid_object( xmax ) ) {
      PLA_Warning( "Invalid object xmax" );
      value--;
    }

   PLA_Obj_objtype( xmax, &objtype );
    if ( objtype != PLA_MSCALAR ){
      PLA_Warning( "Invalid objtype for xmax" );
      value--;
    }      

    PLA_Obj_global_length( xmax, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid global length for xmax" );
      value--;
    }      

    PLA_Obj_global_width( xmax, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid global width for xmax" );
      value--;
    }      

    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    PLA_Obj_objtype( x, &objtype );
    if ( objtype != PLA_MVECTOR )
      PLA_Warning( "Invalid objtype for x" );
      value--;
    }      

  PLA_Obj_global_length( x, &length_x );
  PLA_Obj_global_width( x, &width_x );

  if ( width_x != 1 ){
    PLA_Warning( "x is not of width 1" );
    value--;
  }      
      
  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              x, PLA_ALL_ROWS, PLA_ALL_COLS, &x_cpy );
    PLA_Mscalar_create_conf_to( 
              k, PLA_ALL_ROWS, PLA_ALL_COLS, &k_cpy );
    PLA_Mscalar_create_conf_to( 
              xmax, PLA_ALL_ROWS, PLA_ALL_COLS, &xmax_cpy );

    PLA_Copy( x, x_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Iamax_exit( PLA_Obj x, PLA_Obj k, PLA_Obj xmax )
{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_x, max_y, max_all, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];
  MPI_Datatype
    datatype;

  PLA_Routine_stack_push( "PLA_Iamax_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    int 
      *k_value_p, *k_cpy_value_p;

    PLA_Local_iamax( x_cpy, k_cpy, xmax_cpy );
    PLA_Obj_local_buffer( k_cpy, (void **) &k_cpy_value_p );
    PLA_Obj_local_buffer( k,     (void **) &k_value_p );

    if ( *k_cpy_value_p != *k_value_p ){
      PLA_Warning( "inconsistent k" );
    }

    PLA_Obj_datatype( x, &datatype );
    if ( datatype == MPI_DOUBLE ){
      double 
	*xmax_value_p, *xmax_cpy_value_p;

      PLA_Obj_local_buffer( xmax_cpy, (void **) &xmax_cpy_value_p );
      PLA_Obj_local_buffer( xmax,     (void **) &xmax_value_p );
      
      if ( *xmax_cpy_value_p != *xmax_value_p )
	PLA_Warning( "inconsistent xmax" );
    }
    else if ( datatype == MPI_FLOAT ){
      float
	*xmax_value_p, *xmax_cpy_value_p;

      PLA_Obj_local_buffer( xmax_cpy, (void **) &xmax_cpy_value_p );
      PLA_Obj_local_buffer( xmax,     (void **) &xmax_value_p );
      
      if ( *xmax_cpy_value_p != *xmax_value_p )
	PLA_Warning( "inconsistent xmax" );
    }
    else if ( datatype == MPI_DOUBLE_COMPLEX ){
      double
	*xmax_value_p, *xmax_cpy_value_p;

      PLA_Obj_local_buffer( xmax_cpy, (void **) &xmax_cpy_value_p );
      PLA_Obj_local_buffer( xmax,     (void **) &xmax_value_p );
      
      if ( *xmax_cpy_value_p != *xmax_value_p ){
	PLA_Warning( "inconsistent xmax.real" );
      }

      if ( *(xmax_cpy_value_p+1) != *(xmax_value_p+1) ){
	PLA_Warning( "inconsistent xmax.imaginary" );
      }

    }
    else if ( datatype == MPI_COMPLEX ){
      float
	*xmax_value_p, *xmax_cpy_value_p;

      PLA_Obj_local_buffer( xmax_cpy, (void **) &xmax_cpy_value_p );
      PLA_Obj_local_buffer( xmax,     (void **) &xmax_value_p );
      
      if ( *xmax_cpy_value_p != *xmax_value_p )
	PLA_Warning( "inconsistent xmax" );

      if ( *(xmax_cpy_value_p+1) != *(xmax_value_p+1) )
	PLA_Warning( "inconsistent xmax" );
    }

    PLA_Obj_free( &x_cpy );
    PLA_Obj_free( &k_cpy );
    PLA_Obj_free( &xmax_cpy );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Iamax: memory discrepency" );

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Iamax_exit" ) != 0 )
    PLA_Warning( "PLA_Iamax_exit: history stack corrupted" );

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Iamax" ) != 0 )
    PLA_Warning( "PLA_Iamax: history stack corrupted" );

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_IAMAX_TIMING] += MPI_Wtime() - Iamax_start_time;
  }

  return value;
}

