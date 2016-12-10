/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Show( PLA_Obj Obj )
{
  MPI_Datatype
    datatype;

  PLA_Obj_datatype( Obj, &datatype );

  if ( datatype == MPI_DOUBLE || datatype == MPI_DOUBLE_COMPLEX )
    PLA_Local_show( " ", Obj, "%lf ", " " );
  else if ( datatype == MPI_FLOAT || datatype == MPI_COMPLEX )
    PLA_Local_show( " ", Obj, "%f ", " " );
  else if ( datatype == MPI_INT )
    PLA_Local_show( " ", Obj, "%e ", " " );
  else
    PLA_Abort( "PLA_Show: datatype not yet implemented", __LINE__, __FILE__ );

  return PLA_SUCCESS;
}

