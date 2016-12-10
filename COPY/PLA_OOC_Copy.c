/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

extern int report_case;

/*----------------------------------------------------------------------*/

int PLA_OOC_Copy( PLA_Obj Obj_from, PLA_Obj Obj_to )

/************************************************************************
 
  Copy contents of Obj_from to Obj_to, where one or the other
  is on OOC object

*************************************************************************/

{
  int 
    objtype_from, proj_onto, 
    owner_col, owner_row,
    ooc_from, ooc_to;

  PLA_Obj_ooc( Obj_from, &ooc_from );
  PLA_Obj_ooc( Obj_to, &ooc_to );

  if ( !ooc_from && !ooc_to )
    return PLA_Copy( Obj_from, Obj_to );

  PLA_Obj_objtype( Obj_from, &objtype_from );
  
  switch( objtype_from ){
  case PLA_MSCALAR: 
    if ( ooc_from )
      PLA_Copy_from_ooc_msc( Obj_from, Obj_to );
    else
      PLA_Copy_to_ooc_msc( Obj_from, Obj_to );
    break;
  case PLA_MVECTOR: 
    if ( ooc_from )
      PLA_Copy_from_ooc_mv( Obj_from, Obj_to );
    else
      PLA_Copy_to_ooc_mv( Obj_from, Obj_to );
    break;
  case PLA_MATRIX:  
    if ( ooc_from )
      PLA_Copy_from_ooc_matrix( Obj_from, Obj_to );
    else
      PLA_Copy_to_ooc_matrix( Obj_from, Obj_to );
    break;
  default:
    PLA_Abort ( "OOC copy not supported for objtype", __LINE__, __FILE__ );
  }
  
  return PLA_SUCCESS;
}

int PLA_Copy_from_ooc_msc( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
}

int PLA_Copy_from_ooc_matrix( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  int 
    fd, j, 
    length, width, ldim_from, ldim_to, 
    typesize;
  
  MPI_Datatype
    datatype;

  char
    *buffer_from, *buffer_to;

  PLA_Obj_datatype ( Obj_from, &datatype );
  MPI_Type_size( datatype, &typesize );

  PLA_Obj_local_length( Obj_from, &length );
  PLA_Obj_local_width(  Obj_from, &width );

  PLA_Obj_local_fd( Obj_from, &fd );
  PLA_Obj_local_buffer( Obj_from, (void **) &buffer_from );
  PLA_Obj_local_buffer( Obj_to, (void **) &buffer_to );
  PLA_Obj_local_ldim( Obj_from, &ldim_from );
  PLA_Obj_local_ldim( Obj_to, &ldim_to );

  if ( ldim_from == length && ldim_to == length )
    PLA_Read( ( long ) ( length*width*typesize ), fd, buffer_from, 
	           buffer_to );
  else{
    for ( j=0; j<width; j++ ){
      PLA_Read( ( long ) ( length*typesize ), fd, buffer_from,
		     buffer_to );
      buffer_from += ldim_from * typesize;
      buffer_to +=   ldim_to * typesize;
    }
  }
  
  return PLA_SUCCESS;
}

int PLA_Copy_from_ooc_mv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
}

int PLA_Copy_to_ooc_msc( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
}

int PLA_Copy_to_ooc_matrix( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  int 
    fd, j,
    length, width, ldim_from, ldim_to, typesize;
  
  MPI_Datatype
    datatype;

  char
    *buffer_from, *buffer_to;

  PLA_Obj_datatype ( Obj_from, &datatype );
  MPI_Type_size( datatype, &typesize );

  PLA_Obj_local_length( Obj_from, &length );
  PLA_Obj_local_width(  Obj_from, &width );

  PLA_Obj_local_buffer( Obj_from, ( void ** ) &buffer_from );
  PLA_Obj_local_fd( Obj_to, &fd );
  PLA_Obj_local_buffer( Obj_to, ( void ** ) &buffer_to );
  PLA_Obj_local_ldim( Obj_from, &ldim_from );
  PLA_Obj_local_ldim( Obj_to, &ldim_to );

  if ( ldim_from == length && ldim_to == length ){
    PLA_Write( ( long ) ( length*width*typesize ), buffer_from, fd, 
		buffer_to );
  }
  else{
    for ( j=0; j<width; j++ ){
      PLA_Write( ( long ) ( length*typesize ), buffer_from, fd, 
		  buffer_to );
      buffer_from += ldim_from * typesize;
      buffer_to +=   ldim_to * typesize;
    }
  }
  


  return PLA_SUCCESS;
}

int PLA_Copy_to_ooc_mv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
}

