/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

void create_problem( PLA_Obj A )
/* Compute a random  matrix A */
{
  int  
    me,
    m_A, n_A, ldim_A, 
    i, j;
  static int first_time = TRUE;
  MPI_Datatype
    datatype;

  PLA_Obj_datatype( A, &datatype );
  PLA_Obj_local_length( A, &m_A );
  PLA_Obj_local_width( A, &n_A );
  PLA_Obj_local_ldim( A, &ldim_A );

  MPI_Comm_rank( MPI_COMM_WORLD, &me );
  if ( first_time ){
    srand48( me * 1793 );               /* Seed the random number generator */
    first_time = FALSE;
  }

  if ( datatype == MPI_FLOAT ){
    float *buff_A;

    PLA_Obj_local_buffer( A, (void **) &buff_A );

    for ( j=0; j<n_A; j++ )
      for ( i=0; i<m_A; i++ )
	buff_A[ j*ldim_A + i ] = drand48() * 2.0 - 1.0; 
  }
  else   if ( datatype == MPI_DOUBLE ){
    double *buff_A;

    PLA_Obj_local_buffer( A, (void **) &buff_A );

    for ( j=0; j<n_A; j++ )
      for ( i=0; i<m_A; i++ )
	buff_A[ j*ldim_A + i ] = drand48() * 2.0 - 1.0; 
  }
  else if ( datatype == MPI_COMPLEX ){
    float *buff_A;

    PLA_Obj_local_buffer( A, (void **) &buff_A );

    for ( j=0; j<n_A; j++ )
      for ( i=0; i<m_A*2; i++ )
	buff_A[ j*ldim_A*2 + i ] = drand48() * 2.0 - 1.0; 
  }
  else if ( datatype == MPI_DOUBLE_COMPLEX ){
    double *buff_A;

    PLA_Obj_local_buffer( A, (void **) &buff_A );

    for ( j=0; j<n_A; j++ )
      for ( i=0; i<m_A*2; i++ )
	buff_A[ j*ldim_A*2 + i ] = drand48() * 2.0 - 1.0; 
  }
  else
    PLA_Abort( "create_problem: datatype not yet supported", __LINE__, __FILE__ );

  return;
}
