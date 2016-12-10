/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

#define PLA_MAX_ROUTINE_STACK_DEPTH 200
#define PLA_MAX_NAME_LENGTH 35

int PLA_ERROR_CHECKING = TRUE;
int PLA_CHECK_PARAMETERS = TRUE;
int PLA_CHECK_AGAINST_SEQUENTIAL = FALSE;
int PLA_CHECK_AGAINST_R12 = FALSE;

static char
  routine_stack[ PLA_MAX_ROUTINE_STACK_DEPTH ][ PLA_MAX_NAME_LENGTH ];

static int
  top_of_stack = 0;

/***************************************************************************/

int PLA_Routine_stack_initialize( )
     
/*----------------------------------------------------------------------------

Purpose : Initialize stack that keeps track of routines called

----------------------------------------------------------------------------*/
{
  top_of_stack = 0;

  return PLA_SUCCESS;
}

/***************************************************************************/

int PLA_Routine_stack_push(
       char           *routine_name ) 			    
     
/*----------------------------------------------------------------------------

Purpose : Push name of current routine onto routine stack

IN     routine_name      name of routine being pushed

----------------------------------------------------------------------------*/
{
  if ( top_of_stack <= PLA_MAX_ROUTINE_STACK_DEPTH )
    strcpy( &routine_stack[ top_of_stack ][ 0 ], routine_name );

  /*  { 
    int me, i;
    
    for (i=0; i<top_of_stack; i++)
      printf("..");

    MPI_Comm_rank( MPI_COMM_WORLD, &me );
    printf("%d: Entering %s\n", me, routine_name );
    } */

  top_of_stack++;

  return PLA_SUCCESS;
}

/***************************************************************************/

int PLA_Routine_stack_pop(
       char           *routine_name ) 			    
     
/*----------------------------------------------------------------------------

Purpose : Pop name of current routine from routine stack

OUT    routine_name      name of routine being popped

----------------------------------------------------------------------------*/
{
  top_of_stack--;
  if ( top_of_stack <= PLA_MAX_ROUTINE_STACK_DEPTH )
    strcpy( routine_name, &routine_stack[ top_of_stack ][ 0 ] );
  else
    strcpy( "ROUTINE_STACK_EXCEEDED", routine_name );
  
  /*  { 
    int me, i;
    
    for (i=0; i<top_of_stack; i++)
      printf("..");

    MPI_Comm_rank( MPI_COMM_WORLD, &me );
    printf("%d: Exiting %s\n", me, routine_name );
    } */

  return PLA_SUCCESS;
}


/***************************************************************************/

int PLA_Abort( char *message, int line, char *file )
     
/*----------------------------------------------------------------------------

Purpose : Abort PLAPACK code

IN     message         nerror message

----------------------------------------------------------------------------*/
{
  int me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  printf("aborting!!!\n");
  printf("%3d: PLAPACK fatal error %s.\n in %s line %d\n", 
	 me, message, file, line );

  PLA_Print_stack_history( );
  
  exit( 0 );
}


/***************************************************************************/

int PLA_Warning( char *message )
     
/*----------------------------------------------------------------------------

Purpose : Report PLAPACK warning

IN     message         error message

----------------------------------------------------------------------------*/
{
  int me, i;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  printf("%3d: PLAPACK warning %s.\n", me, message );

  PLA_Print_stack_history( );

  return PLA_SUCCESS;
}

/***************************************************************************/

int PLA_Print_stack_history( )
     
/*----------------------------------------------------------------------------

Purpose : Print routine stack

----------------------------------------------------------------------------*/
{
  int me, i;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  printf("%3d: Stack history:\n", me);

  for ( i=top_of_stack-1; i>=0; i-- ){
    if ( i<=PLA_MAX_ROUTINE_STACK_DEPTH )
      printf( "%3d: %s\n", me, &routine_stack[ i ][ 0 ]);
    else
      printf( "%3d: ...\n", me);
  }
  
  return PLA_SUCCESS;
}

/***************************************************************************/

int PLA_Valid_template( PLA_Template templ )
     
/*----------------------------------------------------------------------------

Purpose : Check if valid template

IN     templ        template to be checked

----------------------------------------------------------------------------*/
{
  return TRUE;
}

/***************************************************************************/

int PLA_Valid_objtype( int objtype )
     
/*----------------------------------------------------------------------------

Purpose : Check if valid object type

IN     objtype      object type to be checked

----------------------------------------------------------------------------*/
{
  int value;

  value = (/* objtype == PLA_VECTOR   || */ objtype == PLA_MVECTOR ||
            objtype == PLA_MATRIX       || objtype == PLA_MSCALAR ||
           /* objtype == PLA_PVECTOR  || */ objtype == PLA_PMVECTOR );

  return value;
}

/***************************************************************************/

int PLA_Valid_object( PLA_Obj obj )
     
/*----------------------------------------------------------------------------

Purpose : Check if valid object

IN     obj        object to be checked

----------------------------------------------------------------------------*/
{
  return PLA_Valid_objtype( obj->objtype );
}

/***************************************************************************/

int PLA_Valid_datatype( MPI_Datatype datatype )
     
/*----------------------------------------------------------------------------

Purpose : Check if valid datatype

IN     datatype      datatype to be checked

----------------------------------------------------------------------------*/
{
  int value;

  value = ( datatype == MPI_INT || datatype == MPI_FLOAT ||
            datatype == MPI_DOUBLE || datatype == MPI_COMPLEX ||
            datatype == MPI_DOUBLE_COMPLEX );

  return value;
}

/***************************************************************************/

int PLA_Compare_local_arrays( int m, int n, char *buffer1, int ldim1,
                                             char *buffer2, int ldim2 )
     
/*----------------------------------------------------------------------------

Purpose : Check if contents of two 2D arrays are identical

IN     m          row dimension
IN     n          column dimension
IN     buffer1    first buffer
IN     ldim1      leading dimenson of first buffer
IN     buffer2    second buffer
IN     ldim2      leading dimenson of second buffer

----------------------------------------------------------------------------*/
{
  int 
    return_value = TRUE,
    i, j;
  char *tempp1, *tempp2;

  for ( j=0; j<n; j++){
    tempp1 = buffer1 + j*ldim1;
    tempp2 = buffer2 + j*ldim2;
    for ( i=0; i<m; i++ ){
      if ( *tempp1 != *tempp2 ) {
	return_value = FALSE;
	break;
      }
    }
    if ( return_value == FALSE ) break;
  }

  return return_value;
}


PLA_Print_objtype( PLA_Obj obj )
{
  int
    objtype;

  PLA_Obj_objtype( obj, &objtype );
  
  switch( objtype ){
  case PLA_MSCALAR:
    printf( "MSCALAR" );
    break;
  case PLA_MATRIX:
    printf( "MATRIX" );
    break;
  case PLA_MVECTOR:
    printf( "MVECTOR");
    break;
  case PLA_PMVECTOR:
    printf( "PMVECTOR");
    break;
  }
}


/***************************************************************************/

int PLA_Valid_side_parameter( int side )
     
/*----------------------------------------------------------------------------

Purpose : Check if valid side parameter

IN     side          side parameter to be checked

----------------------------------------------------------------------------*/
{
  int
    value;

  value = ( side == PLA_SIDE_LEFT || side == PLA_SIDE_RIGHT );

  return value;
}

/***************************************************************************/

int PLA_Valid_trans_parameter( int trans )
     
/*----------------------------------------------------------------------------

Purpose : Check if valid transpose parameter

IN     trans         transpose parameter to be checked

----------------------------------------------------------------------------*/
{
  int
    value;

  value = ( trans == PLA_NO_TRANS ||
	    trans == PLA_TRANS ||
	    trans == PLA_CONJ_TRANS ||
	    trans == PLA_CONJ );

  return value;
}

/***************************************************************************/

int PLA_Valid_uplo_parameter( int uplo )
     
/*----------------------------------------------------------------------------

Purpose : Check if valid uplo parameter

IN     uplo         transpose parameter to be checked

----------------------------------------------------------------------------*/
{
  int
    value;

  value = ( uplo == PLA_LOWER_TRIANGULAR || uplo == PLA_UPPER_TRIANGULAR );

  return value;
}

/***************************************************************************/

int PLA_Valid_diag_parameter( int diag )
     
/*----------------------------------------------------------------------------

Purpose : Check if valid uplo parameter

IN     diag         parameter to be checked

----------------------------------------------------------------------------*/
{
  int
    value;

  value = ( diag == PLA_NONUNIT_DIAG || diag == PLA_UNIT_DIAG );

  return value;
}


/***************************************************************************/

int PLA_Set_error_checking( 
            int error, int parameters, int sequential, int r12 )
/*----------------------------------------------------------------------------

Purpose : Set error checking flags.

IN     error         turn on error checking
IN     parameter     turn on parameter checking
IN     sequential    turn on checking against sequential implementation
IN     r12           turn on checking against PLAPACK version R1.2

Note: is error checking is turned off, all other checking is turned off.

Note: if a parameter is negative, no change is made for that mode.

----------------------------------------------------------------------------*/
{
  if ( error >= 0 ) 
    PLA_ERROR_CHECKING = ( error > 0 );

  if ( parameters >= 0 ) 
    PLA_CHECK_PARAMETERS = ( parameters > 0 );

  if ( sequential >= 0 ) 
    PLA_CHECK_AGAINST_SEQUENTIAL = ( sequential > 0 );

  if ( r12 >= 0 ) 
    PLA_CHECK_AGAINST_R12 = ( r12 > 0 );

  return PLA_SUCCESS;
}

/***************************************************************************/

int PLA_Get_error_checking( 
            int *error, int *parameters, int *sequential, int *r12 )
/*----------------------------------------------------------------------------

Purpose : Set error checking flags.

OUT    error         return if error checking
OUT    parameter     return if parameter checking
OUT    sequential    return if checking against sequential implementation
OUT    r12           return if checking against PLAPACK version R1.2

Note: is error checking is turned off, all other checking is turned off.

Note: if a parameter is negative, no change is made for that mode.

----------------------------------------------------------------------------*/
{
  *error;

  *parameters = PLA_CHECK_PARAMETERS;

  *sequential = PLA_CHECK_AGAINST_SEQUENTIAL;

  *r12 = PLA_CHECK_AGAINST_R12;

  return PLA_SUCCESS;
}
