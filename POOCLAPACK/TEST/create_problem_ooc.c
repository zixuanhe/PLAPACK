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

void create_problem_ooc( int nb_ooc, PLA_Obj A_ooc )
/* Compute a random ooc matrix A_ooc */
{
  int 
    size;

  MPI_Datatype
    datatype;

  PLA_Obj
    A_in_1 = NULL, A_ooc_1 = NULL, A_ooc_cur = NULL;

  PLA_Obj_view_all( A_ooc, &A_ooc_cur );

  while ( TRUE ){
    PLA_Obj_global_width( A_ooc_cur, &size );
    size = min( size, nb_ooc );
    if ( size == 0 ) break;

    PLA_Obj_vert_split_2( A_ooc_cur, size, &A_ooc_1, &A_ooc_cur );

    PLA_Matrix_create_conf_to( A_ooc_1, &A_in_1 );

    create_problem( A_in_1 );

    PLA_Copy( A_in_1, A_ooc_1 );
  }
}
