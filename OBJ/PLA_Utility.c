/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/***************************************************************************/

int PLA_Show_local_contents   (  PLA_Obj       obj,  char *format )

/*----------------------------------------------------------------------------

Purpose : Show local contents of object

IN     obj               global object
IN     format            format for printing
----------------------------------------------------------------------------*/
{
  int
    me, nprocs, curproc, myrow, mycol,
    m, n, ldim, i_out, i_in, i_in_last, 
                j_out, j_in, j_in_last, 
    nb=5;
  MPI_Datatype 
    datatype;
  MPI_Comm 
    comm = MPI_COMM_NULL;
  PLA_Template 
    templ = NULL;

  PLA_Obj_template( obj, &templ );
  PLA_Temp_comm_all     ( templ, &comm );
  PLA_Temp_comm_all_rank( templ, &me );
  PLA_Temp_comm_all_size( templ, &nprocs );
  PLA_Temp_comm_col_rank( templ, &myrow );
  PLA_Temp_comm_row_rank( templ, &mycol );

  for ( curproc=0; curproc < nprocs; curproc++ ){
    MPI_Barrier( comm );
    if ( me == curproc ){
      PLA_Obj_local_length( obj, &m );
      PLA_Obj_local_width ( obj, &n );
      PLA_Obj_local_ldim  ( obj, &ldim );    
      PLA_Obj_datatype    ( obj, &datatype );

      if ( m * n == 0 )
	printf( "NODE %d: no local entries\n" );
      else {
	if ( datatype == MPI_DOUBLE )
	  {
	    double *buf_obj, *tempp;
	    
	    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

	    for ( i_out=0; i_out<m; i_out+=nb ){
	      i_in_last = min( i_out+nb, m );
	      for (j_out=0; j_out<n; j_out+=nb ){
		j_in_last = min( j_out+nb, n );
		printf( "NODE %d  (%d,%d): ROWS %d-%d COLUMNS %d-%d\n", 
		       curproc, myrow, mycol, 
		       i_out, i_in_last-1, j_out, j_in_last-1 );
		for ( i_in=i_out; i_in<i_in_last; i_in++ ){
		  tempp = buf_obj + j_out*ldim + i_in;
		  for ( j_in=j_out; j_in<j_in_last; j_in++){
		    printf(format, *tempp );
		    tempp += ldim;
		  }
		  printf("\n");
		}
		printf("\n");
	      }
	    }
	    break;
	  }
	else{
	  PLA_Warning( "Datatype not yet implemented" );
	}
      }
    }
  }

  MPI_Barrier( comm );

  return PLA_SUCCESS;
}

