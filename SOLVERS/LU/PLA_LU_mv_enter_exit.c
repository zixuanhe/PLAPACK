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

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_LU_mv_enter( PLA_Obj A_mv, PLA_Obj A_11_msc, PLA_Obj pivots )
{
  int 
    value = PLA_SUCCESS,
    length_A_mv, width_A_mv, length_A_11_msc, 
    width_A_11_msc, length_pivots, width_pivots,
    objtype;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_LU_mv";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_LU_mv_enter" );

  old_size_malloced = PLA_Total_size_malloced( );

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if A_mv is valid mvector */

    if ( A_mv == NULL || !PLA_Valid_object( A_mv ) ) {
      PLA_Warning( "Invalid object A_mv" );
      value--;
    }

    PLA_Obj_objtype( A_mv, &objtype );
    if ( objtype != PLA_MVECTOR ){
      PLA_Warning( "Invalid objtype for A_mv" );
      value--;
    }      

    /* Check if datatype is currently supported */
    PLA_Obj_datatype( A_mv, &datatype );
    if ( datatype != MPI_FLOAT && datatype != MPI_DOUBLE &&
	 datatype != MPI_COMPLEX && datatype != MPI_DOUBLE_COMPLEX ){
      PLA_Warning( "datatype not yet supported" );
      value--;
    }

    PLA_Obj_global_length( A_mv, &length_A_mv );
    PLA_Obj_global_width( A_mv, &width_A_mv );

    /* Check if A_11_msc is valid multiscalar */

    if ( A_11_msc == NULL || !PLA_Valid_object( A_11_msc ) ) {
      PLA_Warning( "Invalid object A_11_msc" );
      value--;
    }

    PLA_Obj_objtype( A_11_msc, &objtype );
    if ( objtype != PLA_MSCALAR ){
      PLA_Warning( "Invalid objtype for A_11_msc" );
      value--;
    }      

    PLA_Obj_global_length( A_11_msc, &length_A_11_msc );
    PLA_Obj_global_width( A_11_msc, &width_A_11_msc );

    /* Check if pivots is valid mscalar */

    if ( pivots == NULL || !PLA_Valid_object( pivots ) ) {
      PLA_Warning( "Invalid object pivots" );
      value--;
    }

    PLA_Obj_objtype( pivots, &objtype );
    if ( objtype != PLA_MSCALAR ){
      PLA_Warning( "Invalid objtype for pivots" );
      value--;
    }      

    /* Check if datatype is currently supported */
    PLA_Obj_datatype( pivots, &datatype );
    if ( datatype != MPI_INT ) {
      PLA_Warning( "pivots must be of type MPI_INT" );
      value--;
    }

    PLA_Obj_global_length( pivots, &length_pivots );
    PLA_Obj_global_width( pivots, &width_pivots );

    if ( length_A_11_msc > length_A_mv ){
      PLA_Warning( "length of A_11_msc greater than that of A_mv" );
      value--;
    }

    if ( width_A_11_msc != width_A_mv ){
      PLA_Warning( "width of A_11_msc not equal to that of A_mv" );
      value--;
    }

    if ( width_A_mv != length_pivots ){
      PLA_Warning( "lenght of pivots not equal to width of A_mv" );
      value--;
    }

    if ( width_pivots != 1 ){
      PLA_Warning( "width of pivots not equal to 1" );
      value--;
    }
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    /*
    PLA_Obj 
      A_tmp = NULL, zero = NULL, one = NULL;
    int 
      length, width;
    PLA_Template
      templ;
    
    PLA_Obj_global_length( A, &length );
    PLA_Obj_global_width( A, &width );
    if ( length != width )
      PLA_Abort( "matrix not square", __LINE__, __FILE__ );

    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_tmp );

    PLA_Obj_template( A, &templ );
    PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
			 length, 1, templ, &y_tmp );
    PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
			 length, 1, templ, &x_tmp );
    create_random_data( x_tmp );

    PLA_Copy( A, A_tmp );

    PLA_Create_constants_conf_to( A, NULL, &zero, &one );
    PLA_Local_gemv( PLA_NO_TRANSPOSE, one, A_tmp, x_tmp, zero, y_tmp );

    PLA_Obj_free( &A_tmp );
    PLA_Obj_free( &zero );
    PLA_Obj_free( &one );
    */
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

int PLA_LU_mv_exit( PLA_Obj A, PLA_Obj pivots )

{
  int value = PLA_SUCCESS,
      size_malloced;
  double 
    max_A, diff,
    PLA_Local_abs_max(), PLA_Local_abs_diff();
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_LU_mv_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
    /*
    PLA_Obj 
      A_tmp = NULL;
    int 
      length, width;
    double
      max_x;
    PLA_Template
      templ;
    
    PLA_Mscalar_create_conf_to( 
              A, PLA_ALL_ROWS, PLA_ALL_COLS, &A_tmp );

    PLA_Obj_global_length( A, &length );
    PLA_Obj_global_width( A, &width );

    PLA_Obj_template( A, &templ );
    PLA_Copy( A, A_tmp );

    PLA_Local_trsv( PLA_LOWER_TRIANGULAR, PLA_NO_TRANSPOSE, 
		     PLA_UNIT_DIAG, A_tmp, y_tmp );

    PLA_Local_trsv( PLA_UPPER_TRIANGULAR, PLA_NO_TRANSPOSE, 
		     PLA_NONUNIT_DIAG, A_tmp, y_tmp );

    max_x = PLA_Local_abs_max( x_tmp );

    diff = PLA_Local_abs_diff( x_tmp, y_tmp );

    if ( diff > 0.000001 * max_x ){
      PLA_Warning( "PLA_LU_mv: large relative error encountered" );
      printf("diff = %le, xmax = %le\n", diff, max_x );
      value--;
    }      

    PLA_Obj_free( &A_tmp );
    PLA_Obj_free( &x_tmp );
    PLA_Obj_free( &y_tmp );
    */
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "PLA_LU_mv: memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


