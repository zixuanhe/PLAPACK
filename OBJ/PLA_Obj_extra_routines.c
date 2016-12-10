/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

#if DEBUG==1
int PLA_Obj_set_local_contents   (
              int      trans,                 int      rows_in_buf,
              int      cols_in_buf,           void     *buf,
              int      leading_dim_buf,       int      stride_buf,
              PLA_Obj obj )

/*----------------------------------------------------------------------------

Purpose : Set the local data from the given object.


IN        trans             indicates whether to transpose data
IN        rows_in_buf       row dimension of data buffer
IN        cols_in_buf       column dimension of data buffer
IN        buf               address of data buffer
IN        leading_dim       leading dimension of data buffer
IN        stride_buf        stride of data buffer where data is put
IN/OUT    obj               global object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    m, n, i, j, k, typesize, ldim, stride;
  char 
    *buf_local, *buf_obj, *tempp_local, *tempp_obj;
  MPI_Datatype 
    datatype;
  
  PLA_Obj_local_length( obj, &m );
  PLA_Obj_local_width ( obj, &n );
  PLA_Obj_local_ldim  ( obj, &ldim );    
  PLA_Obj_datatype    ( obj, &datatype );
  MPI_Type_size( datatype, &typesize );
      
  buf_local = (char *) buf;
  PLA_Obj_local_buffer( obj, (void **) &buf_obj );

  if ( trans == PLA_NO_TRANS ){
    for ( j=0; j<n; j++ ){
      tempp_local = buf_local + j*leading_dim_buf*typesize;
      tempp_obj   = buf_obj   + j*ldim*typesize;
      for ( i=0; i<m*typesize; i++ )
	*tempp_obj++ = *tempp_local++;
    }
  }
  else {
    for ( j=0; j<m; j++ ){
      tempp_local = buf_local + j*leading_dim_buf*typesize;
      tempp_obj   = buf_obj + j*typesize;
      for ( i=0; i<n; i++ ) {
	for ( k=0; k<typesize; k++ )
	  *tempp_obj++ = *tempp_local++;
	tempp_obj += (ldim-1)*typesize;
      }
    }
  }

  return PLA_SUCCESS;
}
#endif

