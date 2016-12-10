/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

static PLA_Obj A_cpy;

int PLA_Tri_red_enter( int uplo, PLA_Obj A, PLA_Obj s, PLA_Obj Q )
{
  int 
    value = PLA_SUCCESS,
    length_A, width_A, length_s, width_s, length_Q, width_Q,
    objtype;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Tri_red";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Tri_red_enter" );

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

    /* Check if matrix dimensions match */
    PLA_Obj_global_length( A, &length_A );
    PLA_Obj_global_width(  A, &width_A );

    if ( length_A != width_A ){
      PLA_Warning( "A is not square" );
      value--;
    }      

    /* Check if datatype is currently supported */
    PLA_Obj_datatype( A, &datatype );
    if ( datatype != MPI_FLOAT && datatype != MPI_DOUBLE ){
      PLA_Warning( "datatype not yet supported" );
      value--;
    }

    /* Check if s is valid vector */

    if ( s == NULL || !PLA_Valid_object( s ) ) {
      PLA_Warning( "Invalid object s" );
      value--;
    }

    PLA_Obj_objtype( s, &objtype );
    if ( objtype != PLA_MVECTOR ){
      PLA_Warning( "Invalid objtype for s" );
      value--;
    }      

    /* Check if matrix dimensions match */
    PLA_Obj_global_length( s, &length_s );
    PLA_Obj_global_width(  s, &width_s );

    if ( length_A != length_s ){
      PLA_Warning( "length of s does not match dimension of A" );
      value--;
    }      

    if ( width_s != 1 ){
      PLA_Warning( "s does not have width 1" );
      value--;
    }      

    /* Check if Q is valid matrix */

    if ( Q == NULL || !PLA_Valid_object( Q ) ) {
      PLA_Warning( "Invalid object Q" );
      value--;
    }

    PLA_Obj_objtype( Q, &objtype );
    if ( objtype != PLA_MATRIX ){
      PLA_Warning( "Invalid objtype for Q" );
      value--;
    }      

    /* Check if matrix dimensions match */
    PLA_Obj_global_length( Q, &length_Q );
    PLA_Obj_global_width(  Q, &width_Q );

    if ( length_Q != width_Q ){
      PLA_Warning( "Q is not square" );
      value--;
    }      

    if ( length_A != length_Q ){
      PLA_Warning( "dimensions of A and Q do not match" );
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

int PLA_Tri_red_exit( int uplo, PLA_Obj A, PLA_Obj s, PLA_Obj Q )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];


  PLA_Routine_stack_push( "PLA_Tri_red_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj 
      A_temp = NULL, Q_temp = NULL, B = NULL, C = NULL,
      minus_one = NULL, zero = NULL, one = NULL;
    MPI_Datatype
      datatype;
    double tol;

    PLA_Obj_datatype( A, &datatype );
    if ( datatype == MPI_DOUBLE || datatype == MPI_DOUBLE_COMPLEX )
      tol = 0.000001;
    else
      tol = 0.0001;

    PLA_Create_constants_conf_to( A, &minus_one, &zero, &one );

    PLA_Obj_set_to_zero_below_first_subdiagonal( A );

    PLA_Mscalar_create_conf_to( A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_temp );
    PLA_Mscalar_create_conf_to( Q, PLA_ALL_ROWS, PLA_ALL_COLS, &Q_temp );
    PLA_Copy( Q, Q_temp );
    PLA_Copy( A, A_temp );

    max_A = PLA_Local_abs_max( A_cpy );

    /* Form Q A Q^T */
    PLA_Mscalar_create_conf_to( A_temp, PLA_ALL_ROWS, PLA_ALL_COLS, &B );
    PLA_Mscalar_create_conf_to( A_temp, PLA_ALL_ROWS, PLA_ALL_COLS, &C );
    /* B = Q A */
    PLA_Local_symm( PLA_SIDE_RIGHT, PLA_LOWER_TRIANGULAR, 
                     one, A_temp, Q_temp, zero, B );
    /* C = A - B Q^T*/
    PLA_Local_copy( A_cpy, C );
    PLA_Local_gemm( PLA_NO_TRANS, PLA_TRANS, minus_one, B, Q_temp, one, C );
    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, C );

    diff = PLA_Local_abs_max( C );

    if ( diff > tol * max_A ){
      PLA_Warning( "PLA_Tri_red: large relative error encountered" );
      value--;
    }      

    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &A_temp );
    PLA_Obj_free( &Q_temp );
    PLA_Obj_free( &B );
    PLA_Obj_free( &C );
    PLA_Obj_free( &minus_one );
    PLA_Obj_free( &zero );
    PLA_Obj_free( &one );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Tri_red: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}



int PLA_Obj_set_to_zero_below_first_subdiagonal( PLA_Obj A )
{
  PLA_Obj
    A21 = NULL;

  PLA_Obj_split_4( A, 1, -1, PLA_DUMMY, PLA_DUMMY,
		              &A21,       PLA_DUMMY );

  PLA_Obj_set_to_zero_below_diagonal( A21 );

  PLA_Obj_free( &A21 );
  
  return PLA_SUCCESS;
}


int PLA_Obj_set_to_zero_below_diagonal( PLA_Obj A )
{
  PLA_Obj
    A_BR = NULL, a_21 = NULL;
  int
    size;

  PLA_Obj_view_all( A, &A_BR );

  while( TRUE ){
    PLA_Obj_global_length( A_BR, &size );
    if ( size == 0 ) break;
    
    PLA_Obj_split_4( A_BR, 1, 1, PLA_DUMMY, PLA_DUMMY,
		                  &a_21,      &A_BR );
    
    PLA_Obj_set_to_zero( a_21 );
  }

  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &a_21 );
  
  return PLA_SUCCESS;
}
