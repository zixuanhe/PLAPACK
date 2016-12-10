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

int PLA_OOC_Gemm( int nb_ooc, int transa, int transb,
		   PLA_Obj alpha, PLA_Obj A_ooc, 
		                   PLA_Obj B_ooc,
		   PLA_Obj beta,  PLA_Obj C )
{
  PLA_Obj 
    A_ooc_L = NULL, A_ooc_1 = NULL, A_in_1 = NULL, one = NULL,
    B_ooc_L = NULL, B_ooc_1 = NULL, B_in_1 = NULL;

  int 
    value = PLA_SUCCESS,
    size;

  PLA_Local_scal( beta, C );   
  PLA_Create_constants_conf_to( C, NULL, NULL, &one );

  PLA_Obj_view_all( A_ooc, &A_ooc_L );
  PLA_Obj_view_all( B_ooc, &B_ooc_L );
  
  while ( TRUE ){
    if ( transa == PLA_NO_TRANS ) 
      PLA_Obj_global_width(  A_ooc_L, &size );
    else                          
      PLA_Obj_global_length( A_ooc_L, &size );

    if ( ( size = min( size, nb_ooc ) ) == 0 ) break;

    if ( transa == PLA_NO_TRANS ) 
      PLA_Obj_vert_split_2( A_ooc_L, size, &A_ooc_1, &A_ooc_L );
    else {
      PLA_Obj_horz_split_2( A_ooc_L, size, &A_ooc_1, 
                                            &A_ooc_L );
    }

    if ( transb == PLA_NO_TRANS ) 
      PLA_Obj_horz_split_2( B_ooc_L, size, &B_ooc_1, 
                                            &B_ooc_L );
    else {
      PLA_Obj_vert_split_2( B_ooc_L, size, &B_ooc_1, &B_ooc_L );
    }

    PLA_Matrix_create_conf_to( A_ooc_1, &A_in_1 );

    PLA_Copy( A_ooc_1, A_in_1 );

    PLA_Matrix_create_conf_to( B_ooc_1, &B_in_1 );
    PLA_Copy( B_ooc_1, B_in_1 );

    PLA_Gemm(transa, transb, alpha, A_in_1, B_in_1, beta, C);

  }

  PLA_Obj_free( &A_ooc_L);   PLA_Obj_free( &A_ooc_1 );
  PLA_Obj_free( &A_in_1 );   PLA_Obj_free( &one );
  PLA_Obj_free( &B_ooc_L);   PLA_Obj_free( &B_ooc_1 );
  PLA_Obj_free( &B_in_1 );  


  return value;
}

