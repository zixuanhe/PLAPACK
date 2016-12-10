/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static PLA_Obj x_cpy = NULL, y_cpy = NULL, alpha_cpy = NULL;

static int old_size_malloced;

int PLA_Dot_enter( PLA_Obj x, PLA_Obj y, PLA_Obj alpha )
{
  int
    value = PLA_SUCCESS,
    size, length_x, width_x, length_y, width_y, 
    align_x, align_y,
    length_alpha, width_alpha,
    objtype;
  char 
    routine_name[ 35 ] = "PLA_Dot";
  MPI_Datatype
    datatype;

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Dot_enter" );

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

    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    PLA_Obj_objtype( x, &objtype );
    if ( objtype != PLA_MVECTOR ){
      PLA_Warning( "Invalid objtype for x" );
      value--;
    }      

    PLA_Obj_global_length( x, &length_x );
    PLA_Obj_global_width( x, &width_x );
    PLA_Obj_global_align_row( x, &align_x );

    if ( width_x != 1 ){
      PLA_Warning( "x is not of width 1" );
      value--;
    }      
      
    /* Check if y is valid vector */

    if ( y == NULL || !PLA_Valid_object( y ) ) {
      PLA_Warning( "Invalid object y" );
      value--;
    }

    PLA_Obj_objtype( y, &objtype );
    if ( objtype != PLA_MVECTOR ){
      PLA_Warning( "Invalid objtype for y" );
      value--;
    }      

    PLA_Obj_global_length( y, &length_y );
    PLA_Obj_global_width( y, &width_y );
    PLA_Obj_global_align_row( y, &align_y );

    if ( width_y != 1 ){
      PLA_Warning( "y is not of width 1" );
      value--;
    }      

    if ( length_x != length_y ){
      PLA_Warning( "x and y not of same length" );
      value--;
    }      

    if ( align_x != align_y ){
      PLA_Warning( "Unaligned mvectors not yet supported" );
      value--;
    }      
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              x, PLA_ALL_ROWS, PLA_ALL_COLS, &x_cpy );
    PLA_Mscalar_create_conf_to( 
              y, PLA_ALL_ROWS, PLA_ALL_COLS, &y_cpy );
    PLA_Mscalar_create_conf_to( 
              alpha, PLA_ALL_ROWS, PLA_ALL_COLS, &alpha_cpy );

    PLA_Copy( x, x_cpy );
    PLA_Copy( y, y_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Dot_exit( PLA_Obj x, PLA_Obj y, PLA_Obj alpha )
{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_x, max_y, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];
  MPI_Datatype
    datatype;

  PLA_Routine_stack_push( "PLA_Dot_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Local_dot( x_cpy, y_cpy, alpha_cpy );

    PLA_Obj_datatype( alpha, &datatype );
    if ( datatype == MPI_DOUBLE ){
      double 
	*alpha_value_p, *alpha_cpy_value_p;

      PLA_Obj_local_buffer( alpha_cpy, (void **) &alpha_cpy_value_p );
      PLA_Obj_local_buffer( alpha,     (void **) &alpha_value_p );
      
      if ( dabs( *alpha_cpy_value_p - *alpha_value_p ) > 
	   0.0000001 * dabs( *alpha_cpy_value_p)  ){
	PLA_Warning( "inconsistent alpha" );
	printf("sequential alpha = %lf parallel alpha = %lf\n",
	       *alpha_cpy_value_p, *alpha_value_p );
      }
    }
    else if ( datatype == MPI_FLOAT ){
      float
	*alpha_value_p, *alpha_cpy_value_p;

      PLA_Obj_local_buffer( alpha_cpy, (void **) &alpha_cpy_value_p );
      PLA_Obj_local_buffer( alpha,     (void **) &alpha_value_p );
      
      if ( dabs( *alpha_cpy_value_p - *alpha_value_p ) > 
	   0.0000001 * dabs( *alpha_cpy_value_p)  )
	PLA_Warning( "inconsistent alpha" );
    }
    else{
      PLA_Warning( "checking of datatype not yet supported" );
    }


    PLA_Obj_free( &x_cpy );
    PLA_Obj_free( &y_cpy );
    PLA_Obj_free( &alpha_cpy );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Dot: memory discrepency" );

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Dot_exit" ) != 0 )
    PLA_Warning( "PLA_Dot_exit: history stack corrupted" );

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Dot" ) != 0 )
    PLA_Warning( "PLA_Dot: history stack corrupted" );

  return value;
}

