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

int PLA_Bi_red_enter( int uplo, PLA_Obj A, PLA_Obj sL, PLA_Obj sR, 
		                PLA_Obj U, PLA_Obj V )
{
  int 
    value = PLA_SUCCESS,
    length_A, width_A, 
    length_sL, width_sL, 
    length_sR, width_sR,     
    length_U, width_U,
    length_V, width_V, objtype;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Bi_red";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Bi_red_enter" );

  old_size_malloced = PLA_Total_size_malloced( );

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if uplo parameter is valid */
    if ( !PLA_Valid_uplo_parameter( uplo ) ){
      PLA_Warning( "Invalid parameter uplo" );
      value--;
    }

    if ( uplo != PLA_UPPER_TRIANGULAR ){
      PLA_Abort( "Only reduction to upper bidiagonal currently supported",
		  __LINE__, __FILE__ );
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

    /* Check if matrix dimensions match */
    PLA_Obj_global_length( A, &length_A );
    PLA_Obj_global_width(  A, &width_A );

    /* Check if datatype is currently supported */
    PLA_Obj_datatype( A, &datatype );
    if ( datatype != MPI_FLOAT && datatype != MPI_DOUBLE ){
      PLA_Warning( "datatype not yet supported" );
      value--;
    }

    /* Check if sL and sR are valid vectors */

    if ( sL == NULL || !PLA_Valid_object( sL ) ) {
      PLA_Warning( "Invalid object sL" );
      value--;
    }

    PLA_Obj_objtype( sL, &objtype );
    if ( objtype != PLA_MVECTOR ){
      PLA_Warning( "Invalid objtype for sL" );
      value--;
    }      

    /* Check if sL and sR are valid vectors */

    if ( sR == NULL || !PLA_Valid_object( sR ) ) {
      PLA_Warning( "Invalid object sR" );
      value--;
    }

    PLA_Obj_objtype( sR, &objtype );
    if ( objtype != PLA_MVECTOR ){
      PLA_Warning( "Invalid objtype for sR" );
      value--;
    }      


    /* Check if matrix dimensions match */
    PLA_Obj_global_length( sL, &length_sL );
    PLA_Obj_global_width(  sL, &width_sL );
    PLA_Obj_global_length( sR, &length_sR );
    PLA_Obj_global_width(  sR, &width_sR );

    if ( length_A < width_A )
      PLA_Abort( "Currently, length of A must be at least as large as width of A",
		  __LINE__, __FILE__ );
      
    if ( length_A != length_sL ){
      PLA_Warning( "length of sL does not match row dimension of A" );
      printf("length_A = %d length_sL = %d\n", length_A, length_sL);
      value--;
    }      

    if ( width_sL != 1 ){
      PLA_Warning( "sL does not have width 1" );
      value--;
    }      

    if ( width_A != length_sR ){
      PLA_Warning( "length of sR does not match col dimension of A" );
      value--;
    }      

    if ( width_sR != 1 ){
      PLA_Warning( "sR does not have width 1" );
      value--;
    }      

    /* Check if U is valid matrix */

    if ( U == NULL || !PLA_Valid_object( U ) ) {
      PLA_Warning( "Invalid object U" );
      value--;
    }

    PLA_Obj_objtype( U, &objtype );
    if ( objtype != PLA_MATRIX ){
      PLA_Warning( "Invalid objtype for U" );
      value--;
    }      

    /* Check if U is valid matrix */

    if ( V == NULL || !PLA_Valid_object( V ) ) {
      PLA_Warning( "Invalid object V" );
      value--;
    }

    PLA_Obj_objtype( V, &objtype );
    if ( objtype != PLA_MATRIX ){
      PLA_Warning( "Invalid objtype for V" );
      value--;
    }      

    /* Check if matrix dimensions match */
    PLA_Obj_global_length( U, &length_U );
    PLA_Obj_global_width(  U, &width_U );
    PLA_Obj_global_length( V, &length_V );
    PLA_Obj_global_width(  V, &width_V );

    if ( length_U != width_U ){
      PLA_Warning( "U is not square" );
      value--;
    }      

    if ( length_A != length_U ){
      PLA_Warning( "dimensions of A and U do not match" );
      value--;
    }      
    if ( length_V != width_V ){
      PLA_Warning( "V is not square" );
      value--;
    }      

    if ( width_A != length_V ){
      PLA_Warning( "dimensions of A and V do not match" );
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

int PLA_Bi_red_exit( int uplo, PLA_Obj A, PLA_Obj sL, PLA_Obj sR, 
		     PLA_Obj U, PLA_Obj V )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];


  PLA_Routine_stack_push( "PLA_Bi_red_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj 
      A_temp = NULL, A_12_temp = NULL, U_temp = NULL, V_temp = NULL,
      B = NULL, C = NULL,
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

    PLA_Mscalar_create_conf_to( A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_temp );
    PLA_Mscalar_create_conf_to( U, PLA_ALL_ROWS, PLA_ALL_COLS, &U_temp );
    PLA_Mscalar_create_conf_to( V, PLA_ALL_ROWS, PLA_ALL_COLS, &V_temp );
    PLA_Copy( U, U_temp );
    PLA_Copy( V, V_temp );
    PLA_Copy( A, A_temp );

    PLA_Set_triang_to_zero( PLA_UPPER_TRIANGULAR, PLA_NONUNIT_DIAG, A_temp );
    /*    PLA_Obj_split_4( A_temp, -1, 1, PLA_DUMMY, &A_12_temp,
	  PLA_DUMMY, PLA_DUMMY ); */
    PLA_Obj_vert_split_2( A_temp, 1, PLA_DUMMY, &A_12_temp );
    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, A_12_temp );

    max_A = PLA_Local_abs_max( A_cpy );

    /* Form Q A_temp Q^T */
    PLA_Mscalar_create_conf_to( A_temp, PLA_ALL_ROWS, PLA_ALL_COLS, &B );
    PLA_Mscalar_create_conf_to( A_temp, PLA_ALL_ROWS, PLA_ALL_COLS, &C );
    /* B = U A */
    PLA_Local_gemm( PLA_NO_TRANSPOSE, PLA_NO_TRANSPOSE,
                     one, U_temp, A_temp, zero, B );
    /* C = A_orig - B V^T - A_orig - U A_temp V^T */
    PLA_Local_copy( A_cpy, C );
    PLA_Local_gemm( PLA_NO_TRANS, PLA_TRANS, minus_one, B, V_temp, one, C );
    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, C );
    PLA_Set_triang_to_zero( PLA_UPPER_TRIANGULAR, PLA_NONUNIT_DIAG, C );

    diff = PLA_Local_abs_max( C );

    if ( diff > tol * max_A ){
      PLA_Warning( "PLA_Bi_red: large relative error encountered" );
      value--;
    }      

    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &A_temp );
    PLA_Obj_free( &U_temp );
    PLA_Obj_free( &V_temp );
    PLA_Obj_free( &B );
    PLA_Obj_free( &C );
    PLA_Obj_free( &minus_one );
    PLA_Obj_free( &zero );
    PLA_Obj_free( &one );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Bi_red: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}
