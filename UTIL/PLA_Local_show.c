/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_show( char *head, PLA_Obj Obj, char *format, char *tail )
{
  MPI_Datatype
    datatype;
  int
    i, j,
    length, width, ldim,
    myrow, mycol, me;
  void 
    *buffer;
  PLA_Template
    templ;

  PLA_Obj_datatype( Obj, &datatype );
  PLA_Obj_local_length( Obj, &length );
  PLA_Obj_local_width ( Obj, &width );
  PLA_Obj_local_buffer( Obj, &buffer );
  PLA_Obj_local_ldim  ( Obj, &ldim );

  PLA_Obj_template( Obj, &templ );
  PLA_Temp_comm_col_rank( templ, &myrow );
  PLA_Temp_comm_row_rank( templ, &mycol );
  PLA_Temp_comm_all_rank( templ, &me );

  printf( "%d %d (%3d): %s\n", myrow, mycol, me, head );

  if ( MPI_FLOAT == datatype ){
    float *bufp;

    for (i=0; i<length; i++ ){
      bufp = ( (float *) buffer ) + i;
      for (j=0; j<width; j++ ){
	printf( format, *bufp );
	bufp += ldim;
      }
      printf("\n");
    }
  }
  else   if ( MPI_DOUBLE == datatype ){
    double *bufp;

    for (i=0; i<length; i++ ){
      bufp = ( (double *) buffer ) + i;
      for (j=0; j<width; j++ ){
	printf( format, *bufp );
	bufp += ldim;
      }
      printf("\n");
    }
  }
  else   if ( MPI_COMPLEX == datatype ){
    float *bufp;

    for (i=0; i<length; i++ ){
      bufp = ( (float *) buffer ) + i*2;
      for (j=0; j<width; j++ ){
	printf( format, *bufp );
	printf( format, *(bufp+1) );
	bufp += ldim*2;
      }
      printf("\n");
    }
  }
  else   if ( MPI_DOUBLE_COMPLEX == datatype ){
    double *bufp;

    for (i=0; i<length; i++ ){
      bufp = ( (double *) buffer ) + i*2;
      for (j=0; j<width; j++ ){
	printf( format, *bufp );
	printf( format, *(bufp+1) );
	bufp += ldim*2;
      }
      printf("\n");
    }
  }
  else   if ( MPI_INT == datatype ){
    int *bufp;

    for (i=0; i<length; i++ ){
      bufp = ( (int *) buffer ) + i;
      for (j=0; j<width; j++ ){
	printf( format, *bufp );
	bufp += ldim;
      }
      printf("\n");
    }
  }
  else 
    PLA_Warning( "PLA_Local_show: datatype not supported\n" );

  printf("\n");

  printf( "%s\n", tail );

  return PLA_SUCCESS;
}
    

