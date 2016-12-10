/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

static double Local_BLAS1_start_time;

/***************************************************************************/

int PLA_Local_copy_enter( PLA_Obj x, PLA_Obj y )
{
  int 
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    proj_onto;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_copy";

  if ( PLA_TIMINGS ){
    Local_BLAS1_start_time = MPI_Wtime();
  }

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_copy_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){

  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_copy_exit( PLA_Obj x, PLA_Obj y )
{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_copy_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_COPY_TIMING] += MPI_Wtime() - Local_BLAS1_start_time;
  }

  return value;
}

/***************************************************************************/

int PLA_Local_iamax_enter( PLA_Obj x, PLA_Obj k, PLA_Obj xmax )
{
  int 
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    proj_onto;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_iamax";

  if ( PLA_TIMINGS ){
    Local_BLAS1_start_time = MPI_Wtime();
  }

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_iamax_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){

  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_iamax_exit( PLA_Obj x, PLA_Obj k, PLA_Obj xmax )
{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_iamax_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_IAMAX_TIMING] += MPI_Wtime() - Local_BLAS1_start_time;
  }

  return value;
}

/***************************************************************************/

int PLA_Local_scal_enter( PLA_Obj alpha, PLA_Obj x )
{
  int 
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    proj_onto;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_scal";

  if ( PLA_TIMINGS ){
    Local_BLAS1_start_time = MPI_Wtime();
  }

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_scal_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){

  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_scal_exit( PLA_Obj alpha, PLA_Obj x )
{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_scal_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_SCAL_TIMING] += MPI_Wtime() - Local_BLAS1_start_time;
  }

  return value;
}

/***************************************************************************/

int PLA_Local_inv_scal_enter( PLA_Obj alpha, PLA_Obj x )
{
  int 
    value = PLA_SUCCESS,
    size, length_A, width_A, length_x, width_x, length_y, width_y,
    proj_onto;
  
  MPI_Datatype
    datatype;

  char 
    routine_name[ 35 ] = "PLA_Local_inv_scal";

  if ( PLA_TIMINGS ){
    Local_BLAS1_start_time = MPI_Wtime();
  }

  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Local_inv_scal_enter" );

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){

  }

  PLA_Routine_stack_pop( routine_name );
				 
  return value;
}


int PLA_Local_inv_scal_exit( PLA_Obj alpha, PLA_Obj x )
{
  int value = PLA_SUCCESS,
      size_malloced;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Local_inv_scal_exit" );

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced )
    PLA_Warning( "Memory discrepency" );
  
  PLA_Routine_stack_pop( routine_name );

  PLA_Routine_stack_pop( routine_name );

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_LOCAL_SCAL_TIMING] += MPI_Wtime() - Local_BLAS1_start_time;
  }

  return value;
}

