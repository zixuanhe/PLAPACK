/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_equal_zero ( PLA_Obj Obj )
{
  int
    value = TRUE;
  MPI_Datatype
    datatype;
  int
    i, j,
    length, width, ldim;
  void 
    *buffer;

  PLA_Obj_datatype( Obj, &datatype );
  PLA_Obj_local_length( Obj, &length );
  PLA_Obj_local_width ( Obj, &width );
  PLA_Obj_local_buffer( Obj, &buffer );
  PLA_Obj_local_ldim  ( Obj, &ldim );

  if ( MPI_FLOAT == datatype ){
    float *bufp, f_zero = 0.0;

    for (j=0; j<width; j++ ){
      bufp = ( (float *) buffer ) + j*ldim;
      for (i=0; i<length; i++ ){
	if ( *bufp != f_zero ) {
	  value = FALSE;
	  break;
	}
	bufp++;
      }
      if ( !value ) break;
    }
  }
  else   if ( MPI_DOUBLE == datatype ){
    double *bufp, d_zero = 0.0;

    for (j=0; j<width; j++ ){
      bufp = ( (double *) buffer ) + j*ldim;
      for (i=0; i<length; i++ ){
	if ( *bufp != d_zero ) {
	  value = FALSE;
	  break;
	}
	bufp++;
      }
      if ( !value ) break;
    }
  }
  else   if ( MPI_COMPLEX == datatype ){
    float *bufp, d_zero = 0.0;

    for (j=0; j<width; j++ ){
      bufp = ( (float *) buffer ) + j*ldim*2;
      for (i=0; i<length; i++ ){
	if ( *bufp != d_zero || *(bufp+1) != d_zero ) {
	  value = FALSE;
	  break;
	}
	bufp++; bufp++;
      }
      if ( !value ) break;
    }
  }
  else   if ( MPI_DOUBLE_COMPLEX == datatype ){
    double *bufp, d_zero = 0.0;

    for (j=0; j<width; j++ ){
      bufp = ( (double *) buffer ) + j*ldim*2;
      for (i=0; i<length; i++ ){
	if ( *bufp != d_zero || *(bufp+1) != d_zero ) {
	  value = FALSE;
	  break;
	}
	bufp++; bufp++;
      }
      if ( !value ) break;
    }
  }
  else 
    PLA_Warning( "PLA_Local_equal_zero: datatype not supported\n" );

  return value;
}
