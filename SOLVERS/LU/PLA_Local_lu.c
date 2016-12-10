/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_lu( PLA_Obj A, PLA_Obj pivots )
{
  int 
    m, n, ldim;
  void 
    *buf_A, *buf_pivots;
  MPI_Datatype
    datatype;

  PLA_Obj_local_length( A, &m );
  PLA_Obj_local_width( A, &n );
  if ( m == 0 || n == 0 )
    return PLA_SUCCESS;

  PLA_Obj_local_buffer( A, &buf_A );
  PLA_Obj_local_buffer( pivots, &buf_pivots );
  PLA_Obj_local_ldim( A, &ldim );
  PLA_Obj_datatype( A, &datatype );

  if ( datatype == MPI_DOUBLE ){
    PLA_dlu( &m, &n, (double *) buf_A, &ldim, (int *) buf_pivots );
  }
  else if ( datatype == MPI_FLOAT ){
    PLA_slu( &m, &n, (float *) buf_A, &ldim, (int *) buf_pivots );
  }
  else if ( datatype == MPI_DOUBLE_COMPLEX ){
    PLA_zlu( &m, &n, ( PLA_DOUBLE_COMPLEX *) buf_A, 
	      &ldim, (int *) buf_pivots );
  }
  else if ( datatype == MPI_COMPLEX ){
    PLA_clu( &m, &n, ( PLA_COMPLEX *) buf_A, 
	      &ldim, (int *) buf_pivots );
  }
  else{
    PLA_Warning( "Datatype not yet implemented" );
  }

  return PLA_SUCCESS;
}

