/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

static PLA_Obj A_cpy = NULL, s_cpy = NULL;

int PLA_QR_enter( PLA_Obj A, PLA_Obj s )
{
  int 
    value = PLA_SUCCESS,
    length_A, width_A, length_s, width_s,
    objtype;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_QR";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_QR_enter" );

  old_size_malloced = PLA_Total_size_malloced( );

  if ( PLA_CHECK_PARAMETERS ){
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
    if ( datatype != MPI_FLOAT && datatype != MPI_DOUBLE ){
      PLA_Warning( "datatype not yet supported" );
      value--;
    }

    PLA_Obj_global_width( A, &width_A );

    /* Check if s is valid mvector */

    if ( s == NULL || !PLA_Valid_object( s ) ) {
      PLA_Warning( "Invalid object s" );
      value--;
    }

    PLA_Obj_objtype( s, &objtype );
    if ( objtype != PLA_MVECTOR ){
      PLA_Warning( "Invalid objtype for s" );
      value--;
    }      

    /* Check if datatype is currently supported */
    PLA_Obj_datatype( s, &datatype );
    if ( datatype != MPI_FLOAT && datatype != MPI_DOUBLE ){
      PLA_Warning( "datatype not yet supported" );
      value--;
    }

    PLA_Obj_global_length( s, &length_s );
    PLA_Obj_global_width( s, &width_s );

    if ( width_A != length_s ){
      PLA_Warning( "length of s does not match width of A" );
      value--;
    }
    
    if ( width_s != 1 ){
      PLA_Warning( "width of s is not 1" );
      value--;
    }
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_cpy );
    PLA_Mscalar_create_conf_to( 
              s, PLA_ALL_ROWS, PLA_ALL_COLS, &s_cpy );

    PLA_Copy( A, A_cpy );
    PLA_Copy( s, s_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_QR_exit( PLA_Obj A, PLA_Obj s )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];


  PLA_Routine_stack_push( "PLA_QR_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj 
      A_tmp = NULL, s_tmp = NULL;
    double
      max_A;
    
    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_tmp );
    PLA_Mscalar_create_conf_to( 
              s, PLA_ALL_ROWS, PLA_ALL_COLS, &s_tmp );

    PLA_Copy( A, A_tmp );
    PLA_Copy( s, s_tmp );

    /*
    PLA_Local_qr( A_cpy, s_cpy );

    max_A = PLA_Local_abs_max( A_tmp );

    diff = PLA_Local_abs_diff( A_tmp, A_cpy );

    if ( diff > 0.000001 * max_A ){
      PLA_Warning( "PLA_QR: large relative error encountered" );
      printf("diff = %le, max_A = %le\n", diff, max_A );
      value--;
    }      
    */

    PLA_Obj_free( &A_tmp );
    PLA_Obj_free( &s_tmp );
    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &s_cpy );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_QR: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


