/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

static PLA_Obj A_cpy = NULL, B_cpy = NULL;

int PLA_Pos_def_solve_enter( PLA_Obj A, PLA_Obj B )
{
  int 
    value = PLA_SUCCESS,
    length_A, width_A, length_B, width_B,
    objtype;
  
  MPI_Datatype
    datatype_A, datatype_B;

  char 
    routine_name[ 35 ] = "PLA_Pos_def_solve";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Pos_def_solve_enter" );

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
    PLA_Obj_datatype( A, &datatype_A );
    if ( datatype_A != MPI_FLOAT && datatype_A != MPI_DOUBLE 
	 && datatype_A != MPI_COMPLEX && datatype_A != MPI_DOUBLE_COMPLEX ){
      PLA_Warning( "datatype not yet supported" );
      value--;
    }

    /* Check if B is valid object */

    if ( B == NULL || !PLA_Valid_object( B ) ) {
      PLA_Warning( "Invalid object B" );
      value--;
    }

    PLA_Obj_objtype( B, &objtype );
    if ( objtype != PLA_MATRIX && objtype != PLA_MVECTOR ){
      PLA_Warning( "Invalid objtype for B" );
      value--;
    }      

    /* Check if matrix dimensions match */
    PLA_Obj_global_length( B, &length_B );
    PLA_Obj_global_width(  B, &width_B );

    if ( length_B != width_A ){
      PLA_Warning( "width of A does not match lenght of B" );
      value--;
    }   

    if ( ( objtype == PLA_MVECTOR || objtype == PLA_PMVECTOR )
	 && width_B != 1 ){
      PLA_Warning( "B has objtype PLA_Mvector: width must equal 1" );
      value--;
    }

    /* Check if datatypes match */
    PLA_Obj_datatype( B, &datatype_B );
    if ( datatype_A != datatype_B ){
      PLA_Warning( "datatypes don't match" );
      value--;
    }
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_cpy );

    PLA_Mscalar_create_conf_to( 
              B, PLA_ALL_ROWS, PLA_ALL_COLS, &B_cpy );

    PLA_Copy( A, A_cpy );
    PLA_Copy( B, B_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Pos_def_solve_exit( PLA_Obj A, PLA_Obj B )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_B, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];


  PLA_Routine_stack_push( "PLA_Pos_def_solve_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj A_temp = NULL, B_temp = NULL, one = NULL;
    MPI_Datatype
      datatype;
    double tol;
    int trans;

    PLA_Create_constants_conf_to( A, NULL, NULL, &one );

    PLA_Obj_datatype( A, &datatype );
    if ( datatype == MPI_DOUBLE || datatype == MPI_DOUBLE_COMPLEX )
      tol = 0.000001;
    else
      tol = 0.0001;

    PLA_Mscalar_create_conf_to( A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_temp );
    PLA_Mscalar_create_conf_to( B, PLA_ALL_ROWS, PLA_ALL_COLS, &B_temp );
    PLA_Copy( A, A_temp );
    PLA_Copy( B, B_temp );

    max_B = PLA_Local_abs_max( B_cpy );

    PLA_Local_chol( PLA_LOWER_TRIANGULAR, A_cpy );

    if ( datatype == MPI_DOUBLE || datatype == MPI_FLOAT )
      trans = PLA_TRANSPOSE;
    else
      trans = PLA_CONJUGATE_TRANSPOSE;

    PLA_Local_trsm ( PLA_SIDE_LEFT, PLA_LOWER_TRIANGULAR, 
		      PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG,
		      one, A_cpy, B_cpy );

    PLA_Local_trsm ( PLA_SIDE_LEFT, PLA_LOWER_TRIANGULAR, 
		      trans,          PLA_NONUNIT_DIAG,
		      one, A_cpy, B_cpy );

    diff = PLA_Local_abs_diff( B_temp, B_cpy );

    if ( diff > tol * max_B ){
      PLA_Warning( "PLA_Pos_def_solve: large relative error encountered" );
      value--;
    }      

    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &B_cpy );
    PLA_Obj_free( &A_temp );
    PLA_Obj_free( &B_temp );
    PLA_Obj_free( &one );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_Pos_def_solve: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


