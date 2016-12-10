/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/***************************************************************************/

int PLA_Shift   (   PLA_Obj      obj,    MPI_Datatype   datatype,
                     void         *value )

/*----------------------------------------------------------------------------

Purpose : Add value to diagonal of object obj


IN/OUT          obj          linear algebra object
IN              datatype     datatype of value
IN             *value        value to be used for shift

----------------------------------------------------------------------------*/
{
  int
    return_value = PLA_SUCCESS,
    length, width,
    size;

  MPI_Datatype 
    datatype_obj;

  PLA_Obj
    obj_cur = NULL, obj_11 = NULL;

  PLA_Obj_view_all( obj, &obj_cur );
  PLA_Obj_datatype    ( obj, &datatype_obj );

  while ( TRUE ){
    PLA_Obj_global_length( obj_cur, &size );
    if ( size == 0 ) break;

    PLA_Obj_split_4( obj_cur, 1, 1, &obj_11,   PLA_DUMMY,
		                     PLA_DUMMY, &obj_cur );


    PLA_Obj_local_length( obj_11, &length );
    PLA_Obj_local_width( obj_11, &width );
    if ( length == 1 && width == 1 ){
      if ( datatype_obj == MPI_DOUBLE ) {
	double 
	  *buf_obj, d_value;

	PLA_Obj_local_buffer( obj_11, (void **) &buf_obj );

	if ( datatype == MPI_DOUBLE )
	  d_value = *( ( double *) value );
	else if ( datatype == MPI_FLOAT )
	  d_value = ( double ) *( ( float *) value );
	else
	  PLA_Abort( "Datatype not yet implemented", __LINE__, __FILE__ );

	*buf_obj += d_value;
      }
      else if ( datatype_obj == MPI_FLOAT ) {
	float
	  *buf_obj, d_value;

	PLA_Obj_local_buffer( obj_11, (void **) &buf_obj );

	if ( datatype == MPI_DOUBLE )
	  d_value = *( ( double *) value );
	else if ( datatype == MPI_FLOAT )
	  d_value = ( double ) *( ( float *) value );
	else
	  PLA_Abort( "Datatype not yet implemented", __LINE__, __FILE__ );

	*buf_obj += d_value;
      }
      else {
	PLA_Abort( "PLA_Shift: Datatype not yet implemented", __LINE__, __FILE__ );
      }
    }
  }

  PLA_Obj_free( &obj_cur );
  PLA_Obj_free( &obj_11 );

  return return_value;
}

