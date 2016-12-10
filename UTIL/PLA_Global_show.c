/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Global_show( char *head, PLA_Obj Obj, char *format, char *tail )
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
  PLA_Obj
    Obj_cpy = NULL;

  PLA_Mscalar_create_conf_to( Obj, PLA_ALL_ROWS, PLA_ALL_COLS, &Obj_cpy );
  PLA_Copy( Obj, Obj_cpy );
  
  MPI_Comm_rank( MPI_COMM_WORLD, &me );
  if ( me == 0 ){
    PLA_Obj_datatype( Obj_cpy, &datatype );
    PLA_Obj_local_length( Obj_cpy, &length );
    PLA_Obj_local_width ( Obj_cpy, &width );
    PLA_Obj_local_buffer( Obj_cpy, &buffer );
    PLA_Obj_local_ldim  ( Obj_cpy, &ldim );

    PLA_Obj_template( Obj_cpy, &templ );
    PLA_Temp_comm_col_rank( templ, &myrow );
    PLA_Temp_comm_row_rank( templ, &mycol );
    PLA_Temp_comm_all_rank( templ, &me );

    printf( "%s\n", head );

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
      PLA_Warning( "PLA_Global_show: datatype not supported\n" );

    printf("\n");

    printf( "%s\n", tail );
  }

  PLA_Obj_free( &Obj_cpy );
}
    

