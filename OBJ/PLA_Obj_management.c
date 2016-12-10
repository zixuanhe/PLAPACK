/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/***************************************************************************/

int PLA_Obj_attach_buffer ( PLA_Obj obj, void *buffer, int local_ldim, 
                             int user_buffer )

/*----------------------------------------------------------------------------

Purpose : Attach a given data buffer to the object.

IN     obj               object to which buffer is to be attached
IN     buffer            buffer to be used for data
IN     ldim              leading dimension of buffer
IN     user_buffer       indicates whether PLAPACK provided buffer,
                         or user provides buffer.

----------------------------------------------------------------------------*/
{
  if ( obj->base_obj->local_buffer != NULL ){
    printf( "PLA_Obj_attach_buffer: object already has buffer\n");
    exit( 0 );
  }

  obj->base_obj->local_buffer     = buffer;
  obj->base_obj->local_ldim       = local_ldim;
  obj->base_obj->user_buffer      = user_buffer;

  obj->local_buffer               = buffer;

  return PLA_SUCCESS;
}


static pla_base_obj_counter = 0;

PLA_Obj_base pla_get_base_obj ( )
{
  PLA_Obj_base result;
  
  result = (PLA_Obj_base) PLA_malloc( (size_t) sizeof( struct PLA_Obj_base_struct ) );
  
  if ( result == NULL ) {
    printf("Malloc failed in pla_get_base_obj()\n");
    exit( 0 );
  }

  pla_base_obj_counter++;

  return result;
}

int pla_detach_base_obj( PLA_Obj_base obj_base )
{
  if ( --obj_base->number_of_views == 0 ){
    if ( !obj_base->user_buffer ){
      if ( obj_base->local_buffer != NULL )
	PLA_free( obj_base->local_buffer );
    }
    obj_base->local_buffer = NULL;

    PLA_free( obj_base );
    pla_base_obj_counter--;
  }

  return PLA_SUCCESS;
}

static pla_view_counter = 0;

PLA_Obj pla_get_view ( )
{
  PLA_Obj result;
  
  result = (PLA_Obj) PLA_malloc( (size_t) sizeof( struct PLA_Obj_view_struct ) );
  
  if ( result == NULL ) {
    printf("Malloc failed in pla_get_base_obj()\n");
    exit( 0 );
  }

  pla_view_counter++;

  return result;
}

int pla_free_view( PLA_Obj obj )
{
  if ( obj != NULL ) {
    pla_detach_base_obj( obj->base_obj );

    obj->base_obj = NULL;

    pla_view_counter--;

    PLA_free( obj );
  }

  return PLA_SUCCESS; 
}

int pla_ooc_free_view( PLA_Obj obj )
{
  if ( obj != NULL ) {
    pla_ooc_detach_base_obj( obj->base_obj );

    obj->base_obj = NULL;

    pla_view_counter--;

    PLA_free( obj );
  }

  return PLA_SUCCESS; 
}

int pla_ooc_detach_base_obj( PLA_Obj_base obj_base )
{
  if ( --obj_base->number_of_views == 0 ){
    if ( !obj_base->user_buffer ){
      if ( obj_base->local_buffer != NULL )
	PLA_free( obj_base->local_buffer );
    }
    obj_base->local_buffer = NULL;

    PLA_free( obj_base );
    pla_base_obj_counter--;
  }

  return PLA_SUCCESS;
}


/***************************************************************************/

int PLA_Obj_attach_fd ( PLA_Obj obj, int fd, int local_ldim, 
                             int user_buffer )

/*----------------------------------------------------------------------------

Purpose : Attach a given file descriptor to the object.

IN     obj               object to which buffer is to be attached
IN     fd                file descriptor to be attached
IN     ldim              leading dimension of file
IN     user_buffer       indicates whether PLAPACK provided file,
                         or user provides file.

----------------------------------------------------------------------------*/
{
  if ( obj->base_obj->fd != -1 ){
    printf( "PLA_Obj_attach_fd: object already has file descriptor\n");
    exit( 0 );
  }

  obj->base_obj->local_buffer     = 0;
  obj->base_obj->fd     = fd;
  obj->base_obj->local_ldim       = local_ldim;
  obj->base_obj->user_buffer      = user_buffer;

  obj->local_buffer               = 0;

  return PLA_SUCCESS;
}
