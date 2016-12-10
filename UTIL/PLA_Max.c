/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include <math.h>

#define max(x,y) ( (x) > (y) ? (x) : (y) )
#define abs(x)   ( (x) < 0 ? -(x) : (x) )

double PLA_Local_abs_max ( PLA_Obj Obj )
{
  double 
    value = 0.0;
  
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
    float *bufp;

    for (i=0; i<length; i++ ){
      bufp = ( (float *) buffer ) + i;
      for (j=0; j<width; j++ ){
	value = ( double ) max( value, *bufp );
      }
      bufp += ldim;
    }
  }
  else   if ( MPI_DOUBLE == datatype ){
    double *bufp;

    for (j=0; j<width; j++ ){
      bufp = ( (double *) buffer ) + j*ldim;
      for (i=0; i<length; i++ ){
	value = max( value, abs( *bufp ) );
	bufp++;
      }
    }
  }
  else   if ( MPI_COMPLEX == datatype ){
    float *bufp;

    for (j=0; j<width; j++ ){
      bufp = ( (float *) buffer ) + j*ldim*2;
      for (i=0; i<length; i++ ){
	value = max( value, 
               sqrt( (float) (*bufp * *bufp + *( bufp+1 ) * *( bufp+1 ) ) ) );
	bufp++;
	bufp++;
      }
    }
  }
  else   if ( MPI_DOUBLE_COMPLEX == datatype ){
    double *bufp;

    for (j=0; j<width; j++ ){
      bufp = ( (double *) buffer ) + j*ldim*2;
      for (i=0; i<length; i++ ){
	value = max( value, 
		     sqrt( *bufp * *bufp + *( bufp+1 ) * *( bufp+1 ) ) );
	bufp++;
	bufp++;
      }
    }
  }
  else 
    PLA_Warning( "PLA_Local_abs_max: datatype not supported\n" );

  return value;
}
    

double PLA_Local_abs_diff ( PLA_Obj Obj1, PLA_Obj Obj2 )
{
  double 
    value = 0.0;
  
  MPI_Datatype
    datatype;
  int
    i, j,
    length, width, ldim1, ldim2;
  void 
    *buffer1, *buffer2;

  PLA_Obj_datatype( Obj1, &datatype );
  PLA_Obj_local_length( Obj1, &length );
  PLA_Obj_local_width ( Obj1, &width );
  PLA_Obj_local_buffer( Obj1, &buffer1 );
  PLA_Obj_local_buffer( Obj2, &buffer2 );
  PLA_Obj_local_ldim  ( Obj1, &ldim1 );
  PLA_Obj_local_ldim  ( Obj2, &ldim2 ); 

  if ( MPI_FLOAT == datatype ){
    float *buf1p, *buf2p;

    for (i=0; i<length; i++ ){
      buf1p = ( (float *) buffer1 ) + i;
      buf2p = ( (float *) buffer2 ) + i;
      for (j=0; j<width; j++ ){
	value = ( double ) max( value, abs( *buf1p - *buf2p ) );
	buf1p += ldim1;
	buf2p += ldim2;
      }
    }
  }
  else   if ( MPI_DOUBLE == datatype ){
    double *buf1p, *buf2p;

    for (i=0; i<length; i++ ){
      buf1p = ( (double *) buffer1 ) + i;
      buf2p = ( (double *) buffer2 ) + i;
      for (j=0; j<width; j++ ){
	value = ( double ) max( value, abs( *buf1p - *buf2p ) );
	buf1p += ldim1;
	buf2p += ldim2;
      }
    }
  }
  else   if ( MPI_DOUBLE_COMPLEX == datatype ){
    double *buf1p, *buf2p;

    for (j=0; j<width; j++ ){
      buf1p = ( (double *) buffer1 ) + j*ldim1*2;
      buf2p = ( (double *) buffer2 ) + j*ldim2*2;
      for (i=0; i<length; i++ ){
	value = ( double ) max( value, 
	   sqrt( ( *buf1p - *buf2p ) * ( *buf1p - *buf2p ) +
		 ( *( buf1p+1 ) - *( buf2p+1 ) ) * 
                 ( *( buf1p+1 ) - *( buf2p+1 ) ) ) );

	buf1p +=2;
	buf2p +=2;
      }
    }
  }
  else   if ( MPI_COMPLEX == datatype ){
    float *buf1p, *buf2p;

    for (j=0; j<width; j++ ){
      buf1p = ( (float *) buffer1 ) + j*ldim1*2;
      buf2p = ( (float *) buffer2 ) + j*ldim2*2;
      for (i=0; i<length; i++ ){
	value = ( double ) max( value, 
	   sqrt( ( double ) ( ( *buf1p - *buf2p ) * ( *buf1p - *buf2p ) +
		 ( *( buf1p+1 ) - *( buf2p+1 ) ) * 
                 ( *( buf1p+1 ) - *( buf2p+1 ) ) ) ) );

	buf1p +=2;
	buf2p +=2;
      }
    }
  }
  else 
    PLA_Warning( "PLA_Local_abs_diff: datatype not supported\n" );

  return value;
}
    






