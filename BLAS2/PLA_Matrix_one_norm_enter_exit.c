/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

int PLA_Matrix_one_norm_enter( PLA_Obj A, PLA_Obj alpha )
{
  int
    size,
    value = PLA_SUCCESS,
    objtype;
  char 
    routine_name[ 35 ] = "PLA_Matrix_one_norm";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Matrix_one_norm_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){
    /* Check if alpha is valid multiscalar of size 1x1 */

    if ( alpha == NULL || !PLA_Valid_object( alpha ) ) {
      PLA_Warning( "Invalid object alpha" );
      value--;
    }

   PLA_Obj_objtype( alpha, &objtype );
    if ( objtype != PLA_MSCALAR ){
      PLA_Warning( "Invalid objtype for alpha (currently must be MSCALAR)" );
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
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Matrix_one_norm_exit( PLA_Obj A, PLA_Obj alpha )
{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    PLA_Local_abs_max(), PLA_Local_abs_diff(), diff;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Matrix_one_norm_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj
      A_tmp = NULL, alpha_tmp = NULL, alpha_cpy = NULL;

    PLA_Mscalar_create_conf_to( A, PLA_ALL_ROWS, PLA_ALL_COLS, 
				 &A_tmp );
    PLA_Copy( A, A_tmp );

    PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS, 
				 &alpha_tmp );

    PLA_Copy( alpha, alpha_tmp );

    PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS, 
				 &alpha_cpy );

    PLA_Local_matrix_one_norm( A_tmp, alpha_cpy );

    diff = PLA_Local_abs_diff( alpha_tmp, alpha_cpy );

    if ( diff > 0.000001 ){
      PLA_Warning( "PLA_Matrix_norm_one: large absolute error encountered" );
      value--;
    }      

    PLA_Obj_free( &alpha_cpy );
    PLA_Obj_free( &alpha_tmp );
    PLA_Obj_free( &A_tmp );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Matrix_norm_one: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}



int PLA_Local_matrix_one_norm( PLA_Obj A, PLA_Obj alpha )
{
  int 
    length_A, width_A,
    length_alpha, width_alpha,
    j, lda, i_one = 1;

  MPI_Datatype
    datatype;

  PLA_Obj_local_length( A, &length_A );
  PLA_Obj_local_width( A,  &width_A );

  PLA_Obj_local_length( alpha, &length_alpha );
  PLA_Obj_local_width ( alpha, &width_alpha );

  if ( length_alpha == 0 || width_alpha == 0 ) 
    return PLA_SUCCESS;

  if ( length_A == 0 || width_A == 0 ){
    PLA_Obj_set_to_zero( alpha );
    return PLA_SUCCESS;
  }

  PLA_Obj_local_ldim( A, &lda );
  PLA_Obj_datatype ( A, &datatype );

  if ( datatype == MPI_DOUBLE ){
    double
      *buff_A, *buff_alpha,
      col_nrm1;
    int
      j;

    PLA_Obj_local_buffer( A, (void **) &buff_A );
    PLA_Obj_local_buffer( alpha, (void **) &buff_alpha );

    *buff_alpha = 0.0;
    for ( j=0; j<width_A; j++ ){
      col_nrm1 = PLA_dasum( &length_A, buff_A + j*lda, &i_one );
      *buff_alpha = max( *buff_alpha, col_nrm1 );
    }
  }
  else if ( datatype == MPI_FLOAT ){
    float
      *buff_A, *buff_alpha,
      col_nrm1;
    int
      j;

    PLA_Obj_local_buffer( A, (void **) &buff_A );
    PLA_Obj_local_buffer( alpha, (void **) &buff_alpha );

    *buff_alpha = 0.0;
    for ( j=0; j<width_A; j++ ){
      col_nrm1 = PLA_sasum( &length_A, buff_A + j*lda, &i_one );
      *buff_alpha = max( *buff_alpha, col_nrm1 );
    }
  }
  else if ( datatype == MPI_DOUBLE_COMPLEX ){
    PLA_DOUBLE_COMPLEX
      *buff_A;
    double
      *buff_alpha, col_nrm1;
    int
      j;

    PLA_Obj_local_buffer( A, (void **) &buff_A );
    PLA_Obj_local_buffer( alpha, (void **) &buff_alpha );

    *buff_alpha = 0.0;
    for ( j=0; j<width_A; j++ ){
      col_nrm1 = PLA_dzasum( &length_A, buff_A + j*lda, &i_one );
      *buff_alpha = max( *buff_alpha, col_nrm1 );
    }
  }
  else if ( datatype == MPI_COMPLEX ){
    PLA_COMPLEX
      *buff_A;
    float
      *buff_alpha, col_nrm1;
    int
      j;

    PLA_Obj_local_buffer( A, (void **) &buff_A );
    PLA_Obj_local_buffer( alpha, (void **) &buff_alpha );

    *buff_alpha = 0.0;
    for ( j=0; j<width_A; j++ ){
      col_nrm1 = PLA_scasum( &length_A, buff_A + j*lda, &i_one );
      *buff_alpha = max( *buff_alpha, col_nrm1 );
    }
  }
  else
    PLA_Warning( "PLA_Local_matrix_one_norm: datatype not yet implemented" );

  return PLA_SUCCESS;
}
