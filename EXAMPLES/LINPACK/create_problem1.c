/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

void create_problem( PLA_Obj A, PLA_Obj x, PLA_Obj b )
{
  PLA_Obj  zero = NULL,    one = NULL;
  int size, me, nprocs, i, j, fill_blocksize, this_fill_blocksize, type_size;
  double d_one = 1.0, time;

  void *locA;

  double *local_buf;
  int  local_m, local_n, local_ldim, local_stride, global_length;

  MPI_Datatype
    datatype;
  
  PLA_Obj_global_length( A, &size );

  PLA_Obj_datatype ( A, &datatype );
  MPI_Type_size ( datatype, &type_size);

  MPI_Comm_rank( MPI_COMM_WORLD, &me );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  PLA_Create_constants_conf_to( A, NULL, &zero, &one );

  /*  srand48( me * 1793 ); */

  PLA_Obj_local_length( A, &local_m );
  PLA_Obj_local_width(  A, &local_n );
  PLA_Obj_local_buffer( A, (void **) &local_buf );
  PLA_Obj_local_ldim(   A, &local_ldim );

  if ( 0 == me ) 
    printf("Using PLAPACK application interface to create problem\n");

  MPI_Barrier ( MPI_COMM_WORLD);
  time = MPI_Wtime ();

  PLA_API_begin();
  PLA_Obj_API_open(A);
  PLA_Obj_API_open(x);

  if ( me == 0 ){ 
    fill_blocksize = 10;

    locA = (void *) PLA_malloc( type_size * size * fill_blocksize  );
  
    for (j=0;j< size; j+=fill_blocksize) {
      this_fill_blocksize = min( fill_blocksize, size - j);
      for (i=0; i < size*this_fill_blocksize; i++)  {   
	/* This loop determines the values to put into matrix */
	if ( MPI_DOUBLE == datatype )
	  ((double *)locA)[i]=drand48() * 2.0 - 1.0;      
	else if ( MPI_FLOAT == datatype )
	  ((float *)locA)[i]=drand48() * 2.0 - 1.0;      
	else if ( MPI_DOUBLE_COMPLEX == datatype ) {
	  ((double *)locA)[2*i]=drand48() * 2.0 - 1.0;      
	  ((double *)locA)[2*i+1]=drand48() * 2.0 - 1.0;      
	}
	else if ( 0 == me )
	  printf("Unhandled datatype in create_problem()\n");
      }
      PLA_API_axpy_matrix_to_global(size, 
				    this_fill_blocksize, 
				    &d_one, 
				    locA, 
				    size, 
				    A, 
				    0, j );
    }

    /* processor zero alone fills the vector */
    for (i=0; i<size; i++)
      if ( MPI_DOUBLE == datatype )
	((double *)locA)[i]=drand48() * 2.0 - 1.0;      
      else if ( MPI_FLOAT == datatype )
	((float *)locA)[i]=drand48() * 2.0 - 1.0;      
      else if ( MPI_DOUBLE_COMPLEX == datatype ) {
	((double *)locA)[2*i]=drand48() * 2.0 - 1.0;      
	((double *)locA)[2*i+1]=drand48() * 2.0 - 1.0;      
      }
      else if ( 0 == me )
	printf("Unhandled datatype in create_problem()\n");

    PLA_API_axpy_vector_to_global( size, &d_one, locA, 1, x, 0);
  
    PLA_free( locA );
  }  

  PLA_Obj_API_close(A);
  PLA_Obj_API_close(x);
  PLA_API_end(); 
  
  MPI_Barrier ( MPI_COMM_WORLD);
  time = MPI_Wtime () - time;

  if ( 0 == me ) {
    printf("time for problem creation: %e seconds\n", time);
  }

  PLA_Gemv( PLA_NO_TRANS, one, A, x, zero, b ); 

  PLA_Obj_free( &zero );         PLA_Obj_free( &one );
}

