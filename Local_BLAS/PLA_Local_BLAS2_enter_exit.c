/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

static double Local_BLAS2_start_time;

/***************************************************************************/

int PLA_Local_gemv_enter ( int transa, PLA_Obj alpha, PLA_Obj A,
			    PLA_Obj x, PLA_Obj beta, PLA_Obj y )
{
  int 
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    proj_onto;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_gemv";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_gemv_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_TIMINGS ){
    Local_BLAS2_start_time = MPI_Wtime();
  }

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if transa parameter is valid */
    if ( !PLA_Valid_trans_parameter( transa ) ){
      PLA_Warning( "Invalid parameter transa" );
      value--;
    }

    /* Check if A is valid matrix */

    if ( A == NULL || !PLA_Valid_object( A ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    /* Check if matrix dimensions match */
    PLA_Obj_local_length( A, &length_A );
    PLA_Obj_local_width(  A, &width_A );

    /* Check if alpha is valid multiscalar of size 1x1 */

    if ( alpha == NULL || !PLA_Valid_object( alpha ) ) {
      PLA_Warning( "Invalid object alpha" );
      value--;
    }

    PLA_Obj_local_length( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local length for alpha" );
      value--;
    }      

    PLA_Obj_local_width( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local width for alpha" );
      value--;
    }      

    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    /*
    PLA_Obj_project_onto( x, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_local_length( x, &length_x );
      PLA_Obj_local_width ( x, &width_x );
    }
    else {
      PLA_Obj_local_length( x, &width_x );
      PLA_Obj_local_width ( x, &length_x );
    }
    */

    PLA_Obj_local_length( x, &length_x );
    PLA_Obj_local_width ( x, &width_x );
    if ( width_x > 1 ) {
      PLA_Obj_local_length( x, &width_x );
      PLA_Obj_local_width ( x, &length_x );
    }

    if ( width_x != 1 && width_x != 0 ){
      PLA_Warning( "x has illegal local width" );
      value--;
    }      
      
    /* Check if y is valid vector */

    if ( y == NULL || !PLA_Valid_object( y ) ) {
      PLA_Warning( "Invalid object y" );
      value--;
    }

    /*
    PLA_Obj_project_onto( y, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_local_length( y, &length_y );
      PLA_Obj_local_width ( y, &width_y );
    }
    else {
      PLA_Obj_local_length( y, &width_y );
      PLA_Obj_local_width ( y, &length_y );
    }
    */

    PLA_Obj_local_length( y, &length_y );
    PLA_Obj_local_width ( y, &width_y );
    if ( width_y > 1 ) {
      PLA_Obj_local_length( y, &width_y );
      PLA_Obj_local_width ( y, &length_y );
    }
    if ( width_y != 1 && width_y != 0 ){
      PLA_Warning( "y has illegal local width" );
      value--;
    }      
      
    /* Check if beta is valid multiscalar of size 1x1 */

    if ( beta == NULL || !PLA_Valid_object( beta ) ) {
      PLA_Warning( "Invalid object beta" );
      value--;
    }

    PLA_Obj_local_length( beta, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local length for beta" );
      value--;
    }      

    PLA_Obj_local_width( beta, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local width for beta" );
      value--;
    }      

    /* Check if matrix dimensions match vector lengths */

    if ( length_A != 0 && width_A != 0 ){
      if ( transa == PLA_NO_TRANS || transa == PLA_CONJ ){
	if ( length_A != length_y ){
	  PLA_Warning( "length of A does not match length of y" );
	  value--;
	}      
	
	if ( width_A != length_x ){
	  PLA_Warning( "width of A does not match length of x" );
	  printf("length_A = %d width_A  = %d length_x = %d width_x = %d\n",
		 length_A, width_A, length_x, width_x );
	  value--;
	}      
      }
      else {
	if ( length_A != length_x ){
	  PLA_Warning( "length of A does not match length of x" );
	  value--;
	}      

	if ( width_A != length_y ){
	  PLA_Warning( "width of A does not match length of y" );
	  value--;
	}      
      }
    }
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_gemv_exit( int transa, 
	             PLA_Obj alpha, PLA_Obj A, PLA_Obj x,
                     PLA_Obj beta,  PLA_Obj y )

{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_gemv_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_GEMV_TIMING] += MPI_Wtime() - Local_BLAS2_start_time;
  }

  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Local_symv_enter ( int uplo, PLA_Obj alpha, PLA_Obj A,
			    PLA_Obj x, PLA_Obj beta, PLA_Obj y )
{
  int 
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    proj_onto;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_gemv";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_symv_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_TIMINGS ){
    Local_BLAS2_start_time = MPI_Wtime();
  }

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if uplo parameter is valid */
    if ( !PLA_Valid_uplo_parameter( uplo ) ){
      PLA_Warning( "Invalid parameter uplo" );
      value--;
    }

    /* Check if A is valid matrix */

    if ( A == NULL || !PLA_Valid_object( A ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    /* Check if matrix dimensions match */
    PLA_Obj_local_length( A, &length_A );
    PLA_Obj_local_width(  A, &width_A );

    /* Check if alpha is valid multiscalar of size 1x1 */

    if ( alpha == NULL || !PLA_Valid_object( alpha ) ) {
      PLA_Warning( "Invalid object alpha" );
      value--;
    }

    PLA_Obj_local_length( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local length for alpha" );
      value--;
    }      

    PLA_Obj_local_width( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local width for alpha" );
      value--;
    }      

    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    PLA_Obj_project_onto( x, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_local_length( x, &length_x );
      PLA_Obj_local_width ( x, &width_x );
    }
    else {
      PLA_Obj_local_length( x, &width_x );
      PLA_Obj_local_width ( x, &length_x );
    }

    if ( width_x != 1 && width_x != 0 ){
      PLA_Warning( "x has illegal width" );
      value--;
    }      
      
    /* Check if y is valid vector */

    if ( y == NULL || !PLA_Valid_object( y ) ) {
      PLA_Warning( "Invalid object y" );
      value--;
    }

    PLA_Obj_project_onto( y, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_local_length( y, &length_y );
      PLA_Obj_local_width ( y, &width_y );
    }
    else {
      PLA_Obj_local_length( y, &width_y );
      PLA_Obj_local_width ( y, &length_y );
    }

    if ( width_x != 1 && width_x != 0 ){
      PLA_Warning( "y has illegal width" );
      value--;
    }      
      
    /* Check if beta is valid multiscalar of size 1x1 */

    if ( beta == NULL || !PLA_Valid_object( beta ) ) {
      PLA_Warning( "Invalid object beta" );
      value--;
    }

    PLA_Obj_local_length( beta, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local length for beta" );
      value--;
    }      

    PLA_Obj_local_width( beta, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local width for beta" );
      value--;
    }      

    /* Check if matrix dimensions match vector lengths */
    if ( length_A != 0 && width_A != 0 ){
      if ( length_A != length_x ){
	PLA_Warning( "length of A does not match length of x" );
	value--;
      }      
	
      if ( width_A != length_y ){
	PLA_Warning( "width of A does not match length of y" );
	value--;
      }      
    }
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_symv_exit( int uplo,
	             PLA_Obj alpha, PLA_Obj A, PLA_Obj x,
                     PLA_Obj beta,  PLA_Obj y )

{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_symv_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_SYMV_TIMING] += MPI_Wtime() - Local_BLAS2_start_time;
  }

  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Local_hemv_enter ( int uplo, PLA_Obj alpha, PLA_Obj A,
			    PLA_Obj x, PLA_Obj beta, PLA_Obj y )
{
  int 
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    proj_onto;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_hemv";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_hemv_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_TIMINGS ){
    Local_BLAS2_start_time = MPI_Wtime();
  }

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if uplo parameter is valid */
    if ( !PLA_Valid_uplo_parameter( uplo ) ){
      PLA_Warning( "Invalid parameter uplo" );
      value--;
    }

    /* Check if A is valid matrix */

    if ( A == NULL || !PLA_Valid_object( A ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    /* Check if matrix dimensions match */
    PLA_Obj_local_length( A, &length_A );
    PLA_Obj_local_width(  A, &width_A );

    /* Check if alpha is valid multiscalar of size 1x1 */

    if ( alpha == NULL || !PLA_Valid_object( alpha ) ) {
      PLA_Warning( "Invalid object alpha" );
      value--;
    }

    PLA_Obj_local_length( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local length for alpha" );
      value--;
    }      

    PLA_Obj_local_width( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local width for alpha" );
      value--;
    }      

    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    PLA_Obj_project_onto( x, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_local_length( x, &length_x );
      PLA_Obj_local_width ( x, &width_x );
    }
    else {
      PLA_Obj_local_length( x, &width_x );
      PLA_Obj_local_width ( x, &length_x );
    }

    if ( width_x != 1 && width_x != 0 ){
      PLA_Warning( "x has illegal width" );
      value--;
    }      
      
    /* Check if y is valid vector */

    if ( y == NULL || !PLA_Valid_object( y ) ) {
      PLA_Warning( "Invalid object y" );
      value--;
    }

    PLA_Obj_project_onto( y, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_local_length( y, &length_y );
      PLA_Obj_local_width ( y, &width_y );
    }
    else {
      PLA_Obj_local_length( y, &width_y );
      PLA_Obj_local_width ( y, &length_y );
    }

    if ( width_x != 1 && width_x != 0 ){
      PLA_Warning( "y has illegal width" );
      value--;
    }      
      
    /* Check if beta is valid multiscalar of size 1x1 */

    if ( beta == NULL || !PLA_Valid_object( beta ) ) {
      PLA_Warning( "Invalid object beta" );
      value--;
    }

    PLA_Obj_local_length( beta, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local length for beta" );
      value--;
    }      

    PLA_Obj_local_width( beta, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local width for beta" );
      value--;
    }      

    /* Check if matrix dimensions match vector lengths */
    if ( length_A != 0 && width_A != 0 ){
      if ( length_A != length_x ){
	PLA_Warning( "length of A does not match length of x" );
	value--;
      }      
	
      if ( width_A != length_y ){
	PLA_Warning( "width of A does not match length of y" );
	value--;
      }      
    }
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_hemv_exit( int uplo,
	             PLA_Obj alpha, PLA_Obj A, PLA_Obj x,
                     PLA_Obj beta,  PLA_Obj y )

{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_hemv_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_HEMV_TIMING] += MPI_Wtime() - Local_BLAS2_start_time;
  }

  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Local_trsv_enter( int uplo, int transa, int diag, 
			   PLA_Obj A, PLA_Obj x ) 
{
  int
    value = PLA_SUCCESS,
    length_A, width_A, length_x, width_x,
    proj_onto;
  char 
    routine_name[ 35 ] = "PLA_Local_trsv";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_trsv_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_TIMINGS ){
    Local_BLAS2_start_time = MPI_Wtime();
  }

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if uplo parameter is valid */
    if ( !PLA_Valid_uplo_parameter( uplo ) ){
      PLA_Warning( "Invalid parameter uplo" );
      value--;
    }

    /* Check if transpose parameters are valid */
    if ( !PLA_Valid_trans_parameter( transa ) ){
      PLA_Warning( "Invalid parameter transa" );
      value--;
    }

    /* Check if diag parameter is valid */
    if ( !PLA_Valid_diag_parameter( diag ) ){
      PLA_Warning( "Invalid parameter diag" );
      value--;
    }

    /* Check if A is valid matrix, and square */

    if ( A == NULL || !PLA_Valid_object( A ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    PLA_Obj_local_length( A, &length_A );
    PLA_Obj_local_width ( A, &width_A );
    
    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    PLA_Obj_local_length( x, &length_x );
    PLA_Obj_local_width( x, &width_x );

    if ( width_x > 1 && length_x > 0 ){
      PLA_Warning( "x has illegal dimension(s)" );
      value--;
    }      
      
    /* Check if dimensions match */

    if ( length_A != 0 && width_A != 0 ){
      if ( length_A != width_A ){
	PLA_Warning( "A is not a square matrix" );
	value--;
      }      

      if ( ( width_x == 1 && length_A != length_x ) ||
	   ( length_x == 1 && length_A != width_x ) ) {
	PLA_Warning( "length of A does not match length of x" );
	value--;
      }     
    }
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Local_trsv_exit( int uplo, int transa, int diag, 
			  PLA_Obj A, PLA_Obj x ) 

{
  int value = PLA_SUCCESS,
      size_malloced;

  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_trsv_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "memory discrepency" );
  
  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_TRSV_TIMING] += MPI_Wtime() - Local_BLAS2_start_time;
  }

  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Local_trmv_enter( int uplo, int transa, int diag, 
			   PLA_Obj A, PLA_Obj x ) 
{
  int
    value = PLA_SUCCESS,
    length_A, width_A, length_x, width_x,
    proj_onto;
  char 
    routine_name[ 35 ] = "PLA_Local_trmv";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_trmv_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_TIMINGS ){
    Local_BLAS2_start_time = MPI_Wtime();
  }

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if uplo parameter is valid */
    if ( !PLA_Valid_uplo_parameter( uplo ) ){
      PLA_Warning( "Invalid parameter uplo" );
      value--;
    }

    /* Check if transpose parameters are valid */
    if ( !PLA_Valid_trans_parameter( transa ) ){
      PLA_Warning( "Invalid parameter transa" );
      value--;
    }

    /* Check if diag parameter is valid */
    if ( !PLA_Valid_diag_parameter( diag ) ){
      PLA_Warning( "Invalid parameter diag" );
      value--;
    }

    /* Check if A is valid matrix, and square */

    if ( A == NULL || !PLA_Valid_object( A ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    PLA_Obj_local_length( A, &length_A );
    PLA_Obj_local_width ( A, &width_A );
    
    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    PLA_Obj_project_onto( x, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_local_length( x, &length_x );
      PLA_Obj_local_width( x, &width_x );
    }
    else {
      PLA_Obj_local_width( x, &length_x );
      PLA_Obj_local_length( x, &width_x );
    }

    if ( width_x != 1 && width_x != 0 ){
      PLA_Warning( "x has illegal width" );
      value--;
    }      
      
    /* Check if dimensions match */

    if ( length_A != 0 && width_A != 0 ){
      if ( length_A != width_A ){
	PLA_Warning( "A is not a square matrix" );
	value--;
      }      

      if ( length_A != length_x ){
	PLA_Warning( "length of A does not match length of x" );
	value--;
      }     
    }
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}

#define max(x,y) ( (x) > (y) ? (x) : (y) )

int PLA_Local_trmv_exit( int uplo, int transa, int diag, 
			  PLA_Obj A, PLA_Obj x ) 

{
  int value = PLA_SUCCESS,
      size_malloced;

  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_trmv_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "memory discrepency" );
  
  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_TRMV_TIMING] += MPI_Wtime() - Local_BLAS2_start_time;
  }

  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Local_ger_enter ( PLA_Obj alpha, PLA_Obj x, PLA_Obj y, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    proj_onto;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_ger";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_ger_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_TIMINGS ){
    Local_BLAS2_start_time = MPI_Wtime();
  }

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if A is valid matrix */

    if ( A == NULL || !PLA_Valid_object( A ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    /* Check if matrix dimensions match */
    PLA_Obj_local_length( A, &length_A );
    PLA_Obj_local_width(  A, &width_A );

    /* Check if alpha is valid multiscalar of size 1x1 */

    if ( alpha == NULL || !PLA_Valid_object( alpha ) ) {
      PLA_Warning( "Invalid object alpha" );
      value--;
    }

    PLA_Obj_local_length( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local length for alpha" );
      value--;
    }      

    PLA_Obj_local_width( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local width for alpha" );
      value--;
    }      

    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    PLA_Obj_local_width ( x, &width_x );
    PLA_Obj_local_length( x, &length_x );
    
    if ( width_x > 1 && length_x > 1 ){
      PLA_Warning( "x has illegal local dimensions" );
      value--;
    }      
      
    /* Check if y is valid vector */

    if ( y == NULL || !PLA_Valid_object( y ) ) {
      PLA_Warning( "Invalid object y" );
      value--;
    }

    PLA_Obj_local_length( y, &length_y );
    PLA_Obj_local_width ( y, &width_y );

    if ( width_y > 1 && length_y > 1 ){
      PLA_Warning( "y has illegal local dimensions" );
      value--;
    }      
      
    /* Check if matrix dimensions match vector lengths */

    if ( length_A != 0 && width_A != 0 ){
      if ( ( width_x == 1 && length_A != length_x ) ||
	   ( length_x == 1 && length_A != width_x ) ){
	PLA_Warning( "length of A does not match length of x" );
	value--;
      }      

      if ( ( width_y == 1 && width_A != length_y ) ||
	   ( length_y == 1 && width_A != width_y ) ){
	PLA_Warning( "length of A does not match length of y" );
	value--;
      }      
    }
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_ger_exit( PLA_Obj alpha, PLA_Obj x, PLA_Obj y, PLA_Obj A )
{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_ger_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_GER_TIMING] += MPI_Wtime() - Local_BLAS2_start_time;
  }

  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Local_syr_enter ( int uplo, PLA_Obj alpha, PLA_Obj x, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x,
    proj_onto;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_syr";

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_syr_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_TIMINGS ){
    Local_BLAS2_start_time = MPI_Wtime();
  }

  if ( PLA_CHECK_PARAMETERS ){
    /* Check if A is valid matrix */

    if ( A == NULL || !PLA_Valid_object( A ) ) {
      PLA_Warning( "Invalid object A" );
      value--;
    }

    /* Check if matrix dimensions match */
    PLA_Obj_local_length( A, &length_A );
    PLA_Obj_local_width(  A, &width_A );

    /* Check if alpha is valid multiscalar of size 1x1 */

    if ( alpha == NULL || !PLA_Valid_object( alpha ) ) {
      PLA_Warning( "Invalid object alpha" );
      value--;
    }

    PLA_Obj_local_length( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local length for alpha" );
      value--;
    }      

    PLA_Obj_local_width( alpha, &size );
    if ( size != 1 ){
      PLA_Warning( "Invalid local width for alpha" );
      value--;
    }      

    /* Check if x is valid vector */

    if ( x == NULL || !PLA_Valid_object( x ) ) {
      PLA_Warning( "Invalid object x" );
      value--;
    }

    PLA_Obj_project_onto( x, &proj_onto );
    if ( proj_onto == PLA_PROJ_ONTO_COL ){
      PLA_Obj_local_length( x, &length_x );
      PLA_Obj_local_width ( x, &width_x );
    }
    else {
      PLA_Obj_local_length( x, &width_x );
      PLA_Obj_local_width ( x, &length_x );
    }

    if ( width_x != 1 && width_x != 0 ){
      PLA_Warning( "x has illegal local width" );
      value--;
    }      
      
    /* Check if matrix dimensions match vector lengths */

    if ( length_A != 0 && width_A != 0 ){
      if ( length_A != width_A ){
	PLA_Warning( "matrix A is not square" );
	value--;
      }      

      if ( length_A != length_x ){
	PLA_Warning( "length of A does not match length of x" );
	value--;
      }      
    }
  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_syr_exit( PLA_Obj alpha, PLA_Obj x, PLA_Obj y, PLA_Obj A )
{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_syr_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_SYR_TIMING] += MPI_Wtime() - Local_BLAS2_start_time;
  }

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Local_syr_exit" ) != 0 )
    PLA_Warning( "PLA_Local_syr_exit: history stack corrupted" );

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Local_syr" ) != 0 )
    PLA_Warning( "PLA_Local_syr: history stack corrupted" );

  return value;
}

