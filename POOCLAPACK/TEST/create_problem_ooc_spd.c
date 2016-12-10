/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/***************************************************************************
  1 Nov 1999 

  POOCLAPACK: Parallel Out-of-Core Linear Algebra Package

  Out-of-Core extension of PLAPACK: Parallel Linear Algebra Package
		  
  Copyright (c) 1999 Robert van de Geijn and 
  The University of Texas at Austin.

  Unlike PLAPACK, this program is NOT being released as free software; 
  you CANNOT redistribute it and/or modify it without explicit written 
  permission from The University of Texas.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

  Written under the direction of: 
 	Robert van de Geijn, Department of Computer Sciences,
	University of Texas at Austin  78712.    
	rvdg@cs.utexas.edu
 
  POOCLAPACK was written by Wesley Reiley and Robert van de Geijn
 			
***********i****************************************************************/


#include "PLA.h"

void create_problem_ooc_spd( int nb_tile, PLA_Obj A_ooc )
/* Compute a random symmetric positive definite matrix A_ooc */
{
  int 
    size, size_in, i, j, n;

  MPI_Datatype
    datatype;

  PLA_Obj
    A_in_1 = NULL, A_ooc_1 = NULL, A_ooc_cur = NULL,
    A_ooc_B = NULL;

  double
    d_n;
  PLA_Obj_view_all( A_ooc, &A_ooc_cur );
  PLA_Obj_global_length( A_ooc, &n );

  d_n = (double) n;

  j = 0;
  while ( TRUE ){
    PLA_Obj_global_width( A_ooc_cur, &size );

    size = min( size, nb_tile );
    if ( size == 0 ) break;

    PLA_Obj_vert_split_2( A_ooc_cur, size, &A_ooc_B, &A_ooc_cur );
    
    i = 0;
    while ( TRUE ){
      PLA_Obj_global_length( A_ooc_B, &size_in );

      size_in = min( size_in, nb_tile );
      if ( size_in == 0 ) break;

      PLA_Obj_horz_split_2( A_ooc_B, size_in, &A_ooc_1, 
                                               &A_ooc_B );
      
      PLA_Matrix_create_conf_to( A_ooc_1, &A_in_1 );

      create_problem( A_in_1 );

      if ( i == j ) /* diagonal block.  Shift to force matrix to be
                       symmetric positive definite */
	PLA_Shift( A_in_1, MPI_DOUBLE, &d_n );

      PLA_Copy( A_in_1, A_ooc_1 );

      i++;
    }
    j++;
  }

  PLA_Obj_free( &A_ooc_cur );
  PLA_Obj_free( &A_ooc_B );
  PLA_Obj_free( &A_ooc_1 );
  PLA_Obj_free( &A_in_1 );

/*  return PLA_SUCCESS;*/
}
