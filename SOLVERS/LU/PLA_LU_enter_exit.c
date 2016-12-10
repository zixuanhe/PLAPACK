/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

static PLA_Obj A_cpy = NULL, pivots_cpy = NULL;

int PLA_LU_enter( PLA_Obj A, PLA_Obj pivots )
{
  int 
    value = PLA_SUCCESS,
    length_A, width_A, length_pivots, width_pivots,
    objtype;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_LU";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_LU_enter" );

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
    if ( datatype != MPI_FLOAT && datatype != MPI_DOUBLE &&
	 datatype != MPI_COMPLEX && datatype != MPI_DOUBLE_COMPLEX ){
      PLA_Warning( "datatype not yet supported" );
      value--;
    }

    PLA_Obj_global_length( A, &length_A );
    PLA_Obj_global_width( A, &width_A );

    /* Check if pivots is valid mvector */

    if ( pivots == NULL || !PLA_Valid_object( pivots ) ) {
      PLA_Warning( "Invalid object pivots" );
      value--;
    }

    PLA_Obj_objtype( pivots, &objtype );
    if ( objtype != PLA_MVECTOR ){
      PLA_Abort( "Invalid objtype for pivots. \nNote: in R2.0 pivots must be declared as a multivector of size n x 1,\nNOT a mscalar of size 1 x n",
		  __LINE__, __FILE__ );
      value--;
    }      

    /* Check if datatype is currently supported */
    PLA_Obj_datatype( pivots, &datatype );
    if ( datatype != MPI_INT ){
      PLA_Warning( "pivots must be of type MPI_INT" );
      value--;
    }

    PLA_Obj_global_length( pivots, &length_pivots );

    if ( length_pivots != min( length_A, width_A ) ){
      PLA_Warning("length of pivots must match min( length of A, width of A )");
      value--;
    }
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_cpy );
    PLA_Mscalar_create_conf_to( 
              pivots, PLA_ALL_ROWS, PLA_ALL_COLS, &pivots_cpy );

    PLA_Copy( A, A_cpy );
    PLA_Copy( pivots, pivots_cpy );
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_LU_exit( PLA_Obj A, PLA_Obj pivots )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];


  PLA_Routine_stack_push( "PLA_LU_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    PLA_Obj 
      A_tmp = NULL, pivots_tmp = NULL;
    double
      max_A, tol;
    MPI_Datatype
      datatype;
    
    PLA_Obj_datatype( A, &datatype );
    if ( datatype == MPI_DOUBLE || datatype == MPI_DOUBLE_COMPLEX )
      tol = 0.0000001;
    else
      tol = 0.0001;

    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_tmp );
    PLA_Mscalar_create_conf_to( 
              pivots, PLA_ALL_ROWS, PLA_ALL_COLS, &pivots_tmp );

    PLA_Copy( A, A_tmp );
    PLA_Copy( pivots, pivots_tmp );

    PLA_Local_lu( A_cpy, pivots_cpy );

    max_A = PLA_Local_abs_max( A_tmp );

    diff = PLA_Local_abs_diff( A_tmp, A_cpy );

    if ( diff > tol * max_A ){
      PLA_Warning( "PLA_LU: large relative error encountered" );
      printf("diff = %le, max_A = %le\n", diff, max_A );
      value--;
    }      

    PLA_Obj_free( &A_tmp );
    PLA_Obj_free( &pivots_tmp );
    PLA_Obj_free( &A_cpy );
    PLA_Obj_free( &pivots_cpy );
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_LU: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


