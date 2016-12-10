/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

void Process_answer( PLA_Obj x )
{
  /* Declarations */
  int size, me, i, type_size;
  double d_one = 1.0, value1, value2;
  void *locx = NULL;
  MPI_Datatype datatype;

  PLA_Obj_datatype ( x, &datatype);
  MPI_Type_size ( datatype, &type_size);

  /* What is my rank in the communicator? */
  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  /* Begin Application Interface Mode */
  PLA_API_begin();
    /* Open object x for read/write mode */
    PLA_Obj_API_open( x );

    if ( me == 0 ) {
      /* Create space for receiving vector */

      PLA_Obj_global_length( x, &size );
      locx =  (void *) PLA_calloc( size, type_size);

      for (i=0; i<size; i++){
	if ( MPI_DOUBLE == datatype ) {
	  ((double *)locx)[i] = 0.0;
	}
	if ( MPI_DOUBLE_COMPLEX == datatype ) {
	  ((PLA_DOUBLE_COMPLEX *)locx)[i].real = 0.0;
	  ((PLA_DOUBLE_COMPLEX *)locx)[i].imaginary = 0.0;
	}
	if ( MPI_FLOAT == datatype ) {
	  ((float *)locx)[i] = 0.0;
	}
      }

      /* Initiate extraction of contents of x to local buffer */
      PLA_API_axpy_global_to_vector( size, &d_one, x, 0, locx, 1);
    }
  
  /* Wait until data actually arrives */
  PLA_Obj_API_sync( x );
  
  if ( me == 0 ) {
    /* Print first and last entry of vector */
    if ( MPI_DOUBLE == datatype ) {
      value1 = (double) ((double *)locx)[0] ;
      value2 = (double) ((double *)locx)[size-1] ;
    }
    if ( MPI_DOUBLE_COMPLEX == datatype ) {
      value1 = ((PLA_DOUBLE_COMPLEX *)locx)[0].real ;
      value2 = ((PLA_DOUBLE_COMPLEX *)locx)[size-1].real ;
    }
    if ( MPI_FLOAT == datatype ) {
      value1 = (double) ((float *)locx)[0] ;
      value2 = (double) ((float *)locx)[size-1] ;
    }
    printf("first computed entry: %f\n", value1 );
    printf("last  computed entry: %f\n", value2 );

    /* Free buffer */
    PLA_free( locx );
  }

    /* Close object x */
    PLA_Obj_API_close(x); 
  /* Exit API mode */
  PLA_API_end();
}
