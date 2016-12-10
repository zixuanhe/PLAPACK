/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static int old_size_malloced;

static double Copy_start_time;

static double Save_local_copy_time;

int PLA_Copy_enter( PLA_Obj from, PLA_Obj to )
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Copy";

  if ( PLA_TIMINGS ){
    Copy_start_time = MPI_Wtime();
	Save_local_copy_time = PLA_TIMINGS[PLA_LOCAL_COPY_TIMING];
  }
				 
  PLA_Routine_stack_push( routine_name );

  PLA_Routine_stack_push( "PLA_Copy_enter" );

/*
  {
    int objtype;

    printf("enter PLA_Copy: ");

    PLA_Obj_objtype( from, &objtype );
    switch ( objtype ){
    case PLA_MSCALAR: 
      printf("from MSCALAR ");
      break;
    case PLA_MVECTOR:
      printf("from MVECTOR ");
      break;
    case PLA_MATRIX:
      printf("from MATRIX ");
      break;
    case PLA_PMVECTOR:
      printf("from PMVECTOR ");
      break;
    }

    PLA_Obj_objtype( to, &objtype );
    switch ( objtype ){
    case PLA_MSCALAR: 
      printf("to MSCALAR ");
      break;
    case PLA_MVECTOR:
      printf("to MVECTOR ");
      break;
    case PLA_MATRIX:
      printf("to MATRIX ");
      break;
    case PLA_PMVECTOR:
      printf("to PMVECTOR ");
      break;
    }
    printf("\n");
  }
*/

  old_size_malloced = PLA_Total_size_malloced( );
  
  if ( PLA_CHECK_PARAMETERS ){
  }

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
  }

  PLA_Routine_stack_pop( routine_name );

  return value;
}

int PLA_Copy_exit( PLA_Obj from, PLA_Obj to )

{
  int 
    value = PLA_SUCCESS,
    size_malloced, PLA_Total_size_malloced(), strcmp();
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_push( "PLA_Copy_exit" );

  if ( PLA_CHECK_AGAINST_SEQUENTIAL ){
  }

  size_malloced = PLA_Total_size_malloced( );

  if ( size_malloced != old_size_malloced ){
    PLA_Warning( "PLA_Copy: memory discrepency" );
    printf( "%d, %d", size_malloced,
		  old_size_malloced );
  }

  PLA_Routine_stack_pop( routine_name );

  if ( strcmp( routine_name, "PLA_Copy_exit" ) != 0 ){
    PLA_Warning( "PLA_Copy_exit: history stack corrupted" ); 
    printf("routine name = %s \n", routine_name );
  }

  PLA_Routine_stack_pop( routine_name ); 

  if ( strcmp( routine_name, "PLA_Copy" ) != 0 ){
    PLA_Warning( "PLA_Copy: history stack corrupted" ); 
    printf("routine name = %s \n", routine_name );
  }

  if ( PLA_TIMINGS ){
    PLA_TIMINGS[PLA_COPY_TIMING] += MPI_Wtime() - Copy_start_time;
	PLA_TIMINGS[PLA_LOCAL_COPY_TIMING] = Save_local_copy_time;
  }

  /*  printf("exiting PLA_Copy\n"); */

  return value;
}



