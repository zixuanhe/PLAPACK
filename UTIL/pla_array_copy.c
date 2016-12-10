/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include "mpi.h"

#define FROM(i,j,lda) ( dfrom[ ((j)-1)*lda + (i)-1 ] )
#define TO(i,j,lda)   ( dto[ ((j)-1)*lda + (i)-1 ] )

void pla_array_copy_new( MPI_Datatype datatype,
                     int m, int n, char *from, int ldim_from, 
                                   char *to,   int ldim_to )
{
  int i, j, length;
  char *from_cur, *to_cur;
  
  MPI_Type_size( datatype, &length);
  m *= length;
  ldim_from *= length;
  ldim_to *= length;

  if ( ldim_from == m && ldim_to == m )
    memcpy( to, from, m*n );
  else{
    from_cur = from;
    to_cur = to;
    for ( i=0; i<n; i++ ){
      memcpy( to_cur, from_cur, m );
      to_cur += ldim_to;
      from_cur += ldim_from;
    }
  }
}
