/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_OOC_Obj_free_enter   (PLA_Obj   *obj)

/*----------------------------------------------------------------------------

Purpose : Linear algebra object destructor.


IN/OUT     object      object to be freed
-----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}


int PLA_OOC_Obj_free_exit   (PLA_Obj   *obj)

/*----------------------------------------------------------------------------

Purpose : Linear algebra object destructor.


IN/OUT     obj      object to be freed
-----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_ooc_enter  ( PLA_Obj     obj,    int *ooc )

/*----------------------------------------------------------------------------

Purpose : Extract object type from object.


IN         obj        object to be queried
OUT        ooc        object type of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_ooc";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( ooc == NULL ){
      PLA_Warning( "ooc has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_ooc_exit  ( PLA_Obj     obj,    int *ooc )

/*----------------------------------------------------------------------------

Purpose : Extract object type from object.


IN         obj        object to be queried
OUT        ooc        object type of object

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Obj_local_fd_enter  ( PLA_Obj     obj,    int *fd )

/*----------------------------------------------------------------------------

Purpose : Extract file descriptor from object.


IN         obj        object to be queried
OUT        fd         file descriptor

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_local_fd";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( fd == NULL ){
      PLA_Warning( "fd has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_local_fd_exit  ( PLA_Obj     obj,    int *fd )

/*----------------------------------------------------------------------------

Purpose : Extract file descriptor from object.


IN         obj        object to be queried
OUT        fd         file descriptor

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_pop( routine_name );

  return value;
}
