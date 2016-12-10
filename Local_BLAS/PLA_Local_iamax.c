/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include <mpi.h>

/****************************************************************************
int PLA_Local_iamax( PLA_Obj x, PLA_Obj k, PLA_Obj xmax)
*****************************************************************************/
int PLA_Local_iamax( PLA_Obj x, PLA_Obj k, PLA_Obj xmax )
{
  int
    value = PLA_SUCCESS,
    local_length_x, local_width_x, size, i_one = 1;
  void 
    *x_buf, *xmax_buf;
  int 
    *k_buf;
  MPI_Datatype
    datatype_x;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_iamax_enter( x, k, xmax );

  PLA_Obj_local_length(x, &local_length_x);
  PLA_Obj_local_width (x, &local_width_x);
  PLA_Obj_datatype( x, &datatype_x );

  PLA_Obj_local_buffer(x,    (void **) &x_buf);
  PLA_Obj_local_buffer(xmax, (void **) &xmax_buf);
  PLA_Obj_local_buffer(k,    (void **) &k_buf);

  size = (local_length_x * local_width_x);

  if( MPI_DOUBLE == datatype_x ){
    if ( size > 0 ){
      *k_buf = PLA_idamax( &local_length_x, (double*) x_buf, &i_one) - 1;
      *((double *) xmax_buf) = *( (double *)x_buf + *k_buf );
    }
    else{
      *k_buf = 0;
      *((double *) xmax_buf ) = 0.0;
    }
  }
  if( MPI_FLOAT == datatype_x ){
    if ( size > 0 ){
      *k_buf = PLA_isamax( &local_length_x, (float *) x_buf, &i_one) - 1;
      *((float *) xmax_buf) = *( (float *)x_buf + *k_buf );
    }
    else{
      *k_buf = 0;
      *((float *) xmax_buf ) = 0.0;
    }
  }
  if( MPI_DOUBLE_COMPLEX == datatype_x ){
    if ( size > 0 ){
      *k_buf = PLA_izamax( &local_length_x, ( PLA_DOUBLE_COMPLEX *) x_buf, &i_one) - 1;
      *((double *) xmax_buf) = *( (double *)x_buf + *k_buf * 2 );
      *((double *) xmax_buf+1) = *( (double *)x_buf + *k_buf * 2 + 1);
    }
    else{
      *k_buf = 0;
      *((double *) xmax_buf )    = 0.0;
      *((double *) xmax_buf + 1) = 0.0;
    }
  }
  if( MPI_COMPLEX == datatype_x ){
    if ( size > 0 ){
      *k_buf = PLA_icamax( &local_length_x, ( PLA_COMPLEX *) x_buf, &i_one) - 1;
      *((float *) xmax_buf) = *( (float *)x_buf + *k_buf * 2 );
      *((float *) xmax_buf+1) = *( (float *)x_buf + *k_buf * 2 + 1);
    }
    else{
      *k_buf = 0;
      *((float *) xmax_buf )    = 0.0;
      *((float *) xmax_buf + 1) = 0.0;
    }
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Local_iamax_exit( x, k, xmax );

  return value;
}
