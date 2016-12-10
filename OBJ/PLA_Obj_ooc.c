/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/***************************************************************************/

int PLA_OOC_Mvector_create_without_fd (
       MPI_Datatype   datatype,        int            global_length,
       int            global_width,    PLA_Template   templ,
       int            global_align,    PLA_Obj       *new_obj)
     
/*----------------------------------------------------------------------------

Purpose : Create distributed multi-vector, without buffer for data.

IN     datatype          datatype of object
IN     global_length     global length of multi-vector
IN     global_width      global width of multi-vector
IN     templ             template for vector and matrix distribution
IN     global_align      alignment to template
OUT    new_obj           object describing created multi-vector

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;

  value = PLA_Mvector_create_without_buffer( datatype, global_length,
				      global_width, templ, global_align,
				      new_obj );

  if ( value == PLA_SUCCESS ){
    (*new_obj)->base_obj->ooc = TRUE;
    (*new_obj)->base_obj->fd = -1;
  }

  return value;
}

/***************************************************************************/

int PLA_OOC_Matrix_create_without_fd   (
     MPI_Datatype   datatype,          int            global_length,
     int            global_width,      PLA_Template   templ,
     int            global_align_row,  int            global_align_col,
     PLA_Obj       *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create distributed matrix without buffer.

IN     datatype          datatype of object
IN     global_length     global length of matrix
IN     global_width      global width of matrix
IN     templ             template for vector and matrix distribution
IN     global_align_row  row alignment to template
IN     global_align_col  column alignment to template
OUT    new_obj           object describing created matrix

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;

  value = PLA_Matrix_create_without_buffer (
     datatype, global_length, global_width, templ, global_align_row,  
     global_align_col, new_obj );

  if ( value == PLA_SUCCESS ){
    ( *new_obj )->base_obj->ooc = TRUE;
    (*new_obj)->base_obj->fd = -1;
  }

  return value;
}

/***************************************************************************/

int PLA_OOC_Mscalar_create_without_fd  (
       MPI_Datatype   datatype,          int            owner_row,
       int            owner_col,         int            length,
       int            width,             PLA_Template   templ,
       PLA_Obj       *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create distributed multiscalar without buffer.


IN     datatype          datatype of object
IN     owner_row         index of row of nodes with contains owning node(s)
IN     owner_col         index of column of nodes with contains owning node(s)
IN     length            length of mscalar
IN     width             width of mscalar
IN     templ             template for vector and matrix distribution
OUT    new_obj           object describing created multiscalar

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;

  value = PLA_Mscalar_create_without_buffer (
     datatype, owner_row, owner_col, length, width, templ, new_obj );

  if ( value == PLA_SUCCESS ){
    ( *new_obj )->base_obj->ooc = TRUE;
    (*new_obj)->base_obj->fd = -1;
  }

  return value;
}

/***************************************************************************/

int PLA_OOC_Obj_free   (PLA_Obj   *obj)

/*----------------------------------------------------------------------------

Purpose : OOC linear algebra object destructor.

IN/OUT     object      object to be freed

-----------------------------------------------------------------------------*/
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_OOC_Obj_free_enter( obj );

  pla_ooc_free_view( *obj );

  *obj = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_OOC_Obj_free_exit( obj );

  return value;
}

/***************************************************************************/

int PLA_Obj_ooc  ( PLA_Obj     obj,    int *ooc )

/*----------------------------------------------------------------------------

Purpose : Extract whether object is OOC.


IN         obj        object to be queried
OUT        ooc        indicates if object is OOC ( TRUE or FALSE )

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_ooc_enter( obj, ooc );

  /* if ( value == PLA_SUCCESS ) */ 
  *ooc = obj->base_obj->ooc;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_ooc_exit( obj, ooc );

  return value;
}

/***************************************************************************/

int PLA_Obj_local_fd  ( PLA_Obj     obj,    int *fd )

/*----------------------------------------------------------------------------

Purpose : Extract file discriptor

IN         obj        object to be queried
OUT        fd         file descriptor

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_fd_enter( obj, fd );

  /* if ( value == PLA_SUCCESS ) */ 
  *fd = obj->base_obj->fd;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_fd_exit( obj, fd );

  return value;
}


/***************************************************************************/

int PLA_Obj_local_status ( PLA_Obj obj, struct ffsw *status)

/*----------------------------------------------------------------------------

Purpose : Extract status of async IO

IN         obj        object to be queried
OUT        status     status of async IO

----------------------------------------------------------------------------*/
{

  int
    value = PLA_SUCCESS;


  status = &obj->base_obj->ooc_status;


  return value;

}
