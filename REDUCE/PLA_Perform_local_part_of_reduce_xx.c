/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include "PLA_Reduce.h"

extern int report_case;

/*----------------------------------------------------------------------*/

int PLA_Perform_local_part_of_reduce_xx( 
                    int shape, PLA_Obj Obj_from, MPI_Op op, 
                    PLA_Obj alpha, PLA_Obj Obj_to )

/************************************************************************
 
  Perform local part of reduce contents of Obj_from to Obj_to
  Expert version: allows shape and adding to multiple of target

*************************************************************************/

{
  if ( op == MPI_SUM ){
    if ( shape == PLA_SHAPE_GENERAL ){
      PLA_Local_scal( alpha, Obj_to );
      PLA_Local_add( Obj_from, Obj_to );
    }
    else
      PLA_Warning( "Shape not yet supported" );
  }
  else
      PLA_Warning( "Operation not yet supported" );
  
  return PLA_SUCCESS;
}
  

