/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/****************************************************************************
int PLA_Local_absolute_max( PLA_Obj x, PLA_Obj xmax)
*****************************************************************************/
int PLA_Local_absolute_max( PLA_Obj x, PLA_Obj xmax )
{
  int  
    local_length_x, local_width_x, size, i_one = 1;
  void 
    *x_buf, *xmax_buf;
  int 
    index;
  MPI_Datatype
    datatype_x;

  PLA_Obj_local_length(x, &local_length_x);
  PLA_Obj_local_width (x, &local_width_x);
  PLA_Obj_datatype( x, &datatype_x );

  PLA_Obj_local_buffer(x,    (void **) &x_buf);
  PLA_Obj_local_buffer(xmax, (void **) &xmax_buf);

  size = (local_length_x * local_width_x);

  if( MPI_DOUBLE == datatype_x ){
    if ( size > 0 ){
      index = PLA_idamax( &local_length_x, (double*) x_buf, &i_one) - 1;
      *((double *) xmax_buf) = *( (double *)x_buf + index );
    }
    else{
      *((double *) xmax_buf ) = 0.0;
    }
  }
  if( MPI_FLOAT == datatype_x ){
    if ( size > 0 ){
      index = PLA_isamax( &local_length_x, (float *) x_buf, &i_one) - 1;
      *((float *) xmax_buf) = *( (float *)x_buf + index );
    }
    else{
      *((float *) xmax_buf ) = 0.0;
    }
  }
  if( MPI_DOUBLE_COMPLEX == datatype_x ){
    if ( size > 0 ){
      index = PLA_izamax( &local_length_x, ( PLA_DOUBLE_COMPLEX *) x_buf, &i_one) - 1;
      *((double *) xmax_buf) = *( (double *)x_buf + index * 2 );
      *((double *) xmax_buf+1) = *( (double *)x_buf + index * 2 + 1);
    }
    else{
      *((double *) xmax_buf )    = 0.0;
      *((double *) xmax_buf + 1) = 0.0;
    }
  }
  if( MPI_COMPLEX == datatype_x ){
    if ( size > 0 ){
      index = PLA_icamax( &local_length_x, ( PLA_COMPLEX *) x_buf, &i_one) - 1;
      *((float *) xmax_buf) = *( (float *)x_buf + index * 2 );
      *((float *) xmax_buf+1) = *( (float *)x_buf + index * 2 + 1);
    }
    else{
      *((float *) xmax_buf )    = 0.0;
      *((float *) xmax_buf + 1) = 0.0;
    }
  }
}
