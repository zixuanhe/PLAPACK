/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/***************************************************************************
  1 Nov 1999 

  POOCLAPACK: Parallel Out-of-Core Linear Algebra Package

  Out-of-Core extension of PLAPACK: Parallel Linear Algebra Package
		  
  Copyright (c) 1999 Robert van de Geijn and 
  The University of Texas at Austin.

  Unlike PLAPACK, this program is NOT being released as free software; 
  you CANNOT redistribute it and/or modify it without explicit written 
  permission from The University of Texas.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

  Written under the direction of: 
 	Robert van de Geijn, Department of Computer Sciences,
	University of Texas at Austin  78712.    
	rvdg@cs.utexas.edu
 
  POOCLAPACK was written by Wesley Reiley and Robert van de Geijn
 			
***********i****************************************************************/


#include "PLA.h"

static PLA_Obj alpha_cpy = NULL, A_cpy = NULL, 
                beta_cpy  = NULL, C_cpy = NULL;

static int old_size_malloced;

int PLA_OOC_Syrk_enter ( int nb_ooc, int uplo, int transa, 
			  PLA_Obj alpha, PLA_Obj A_ooc,
			  PLA_Obj beta,  PLA_Obj C )
{
  int
    value = PLA_SUCCESS,
    size, length_A, width_A, length_C, width_C,
    objtype, ooc;

  PLA_Obj
    A_in = NULL;

  char 
    routine_name[ 35 ] = "PLA_OOC_Syrk";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_OOC_Syrk_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){
    /* Check if transpose parameter is valid */
    if ( !PLA_Valid_trans_parameter( transa ) ){
      PLA_Warning( "Invalid parameter transa" );
      value--;
    }

    /* Check if uplo parameter is valid */
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

    if ( A_ooc == NULL || !PLA_Valid_object( A_ooc ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    PLA_Obj_objtype( A_ooc, &objtype );
    if ( objtype != PLA_MATRIX ){
      PLA_Warning( "Invalid objtype for A_ooc" );
      value--;
    }      

    PLA_Obj_ooc( A_ooc, &ooc );
    if ( ooc != TRUE ){
      PLA_Warning( "A_ooc not an OOC matrix" );
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

    /* Check if C is valid matrix */

    if ( C == NULL || !PLA_Valid_object( C ) ) {
      PLA_Warning( "Invalid object C" );
      value--;
    }

    PLA_Obj_objtype( C, &objtype );
    if ( objtype != PLA_MATRIX ){
      PLA_Warning( "Invalid objtype for C" );
      value--;
    }      

    PLA_Obj_ooc( C, &ooc );
    if ( ooc != FALSE ){
      PLA_Warning( "C cannot be an OOC matrix" );
      value--;
    }      

    /* Check if matrix dimensions match */
    PLA_Obj_global_length( A_ooc, &length_A );
    PLA_Obj_global_width(  A_ooc, &width_A );
    PLA_Obj_global_length( C, &length_C );
    PLA_Obj_global_width(  C, &width_C );

    if ( length_C != width_C ){
      PLA_Warning( "C is not a square matrix" );
      value--;
    }      

    if ( transa == PLA_NO_TRANS || transa == PLA_CONJ ){
      if ( length_A != length_C ){
	PLA_Warning( "length of A does not match length of C" );
	value--;
      }      
    }
    else {
      if ( width_A != length_C ){
	PLA_Warning( "width of A does not match length of C" );
	value--;
      }      
    }
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              alpha, PLA_ALL_ROWS, PLA_ALL_COLS, &alpha_cpy );
    PLA_Mscalar_create_conf_to( 
              A_ooc, PLA_ALL_ROWS, PLA_ALL_COLS, &A_cpy );
    PLA_Mscalar_create_conf_to( 
              beta, PLA_ALL_ROWS, PLA_ALL_COLS, &beta_cpy );
    PLA_Mscalar_create_conf_to( 
              C, PLA_ALL_ROWS, PLA_ALL_COLS, &C_cpy );


    PLA_Copy( alpha, alpha_cpy );
    PLA_Matrix_create_conf_to( A_ooc, &A_in );
    PLA_Copy( A_ooc, A_in );
    PLA_Copy( A_in, A_cpy );
    PLA_Copy( beta, beta_cpy );
    if ( PLA_Local_equal_zero( beta_cpy ) )
      PLA_Obj_set_to_zero( C_cpy );
    else
      PLA_Copy( C, C_cpy );

    PLA_Obj_free( &A_in );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_OOC_Syrk_exit ( int nb_ooc, int uplo, int transa, 
		     PLA_Obj alpha, PLA_Obj A,
		     PLA_Obj beta,  PLA_Obj C )
{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, max_C, max_all, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  PLA_Obj
    one = NULL;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_OOC_Syrk_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj C_temp = NULL;

    PLA_Mscalar_create_conf_to( C, PLA_ALL_ROWS, PLA_ALL_COLS, &C_temp );
    PLA_Copy( C, C_temp );

    max_A = PLA_Local_abs_max( A_cpy );

    PLA_Local_scal( beta_cpy, C_cpy );

    max_C = PLA_Local_abs_max( C_cpy ); 
    max_all = max( max_A, max_C ); 

    PLA_Create_constants_conf_to( A, NULL, NULL, &one );
    PLA_Local_syrk( uplo, transa, alpha_cpy, A_cpy, one,  C_cpy );

    diff = PLA_Local_abs_diff( C_temp, C_cpy );

    if ( diff > 0.000001 * max_all ){
      PLA_Warning( "PLA_OOC_Syrk: large relative error encountered" );
      value--;
    }      

    PLA_Obj_free( &alpha_cpy );
    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &beta_cpy );
    PLA_Obj_free( &C_cpy );

    PLA_Obj_free( &C_temp );
    PLA_Obj_free( &one );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_OOC_Syrk: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


