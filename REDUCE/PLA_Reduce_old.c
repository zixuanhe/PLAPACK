/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Reduce_old( PLA_Obj from, MPI_Op op, PLA_Obj to )
{
  int 
    objtype_from, objtype_to, owner_row, owner_col, owner;

  PLA_Obj_objtype( from, &objtype_from );
  PLA_Obj_objtype( to,   &objtype_to );

  switch( objtype_from ){
  case PLA_MATRIX:
    printf("Obj from: MATRIX\n");
    break;
  case PLA_PMVECTOR:
    printf("Obj from: PMVECTOR\n");
    break;
  case PLA_MVECTOR:
    printf("Obj from: MVECTOR\n");
    break;
  case PLA_MSCALAR:
    printf("Obj from: MSCALAR\n");

    PLA_Obj_owner_row ( from, &owner );
    if ( owner == PLA_ALL_ROWS ) 
      printf("owner_row = PLA_ALL_ROWS\n");
    else
      printf("owner_row = %d\n", owner);

    PLA_Obj_owner_col ( from, &owner );
    if ( owner == PLA_ALL_COLS ) 
      printf("owner_col = PLA_ALL_COLS\n");
    else
      printf("owner_col = %d\n", owner);

    break;
  default:
    printf("Unknown objtype_from\n");
  }

  switch( objtype_to ){
  case PLA_MATRIX:
    printf("Obj to: MATRIX\n");
    break;
  case PLA_PMVECTOR:
    printf("Obj to: PMVECTOR\n");
    break;
  case PLA_MVECTOR:
    printf("Obj to: MVECTOR\n");
    break;
  case PLA_MSCALAR:
    printf("Obj to: MSCALAR\n");

    PLA_Obj_owner_row ( to, &owner );
    if ( owner == PLA_ALL_ROWS ) 
      printf("owner_row = PLA_ALL_ROWS\n");
    else
      printf("owner_row = %d\n", owner);

    PLA_Obj_owner_col ( to, &owner );
    if ( owner == PLA_ALL_COLS ) 
      printf("owner_col = PLA_ALL_COLS\n");
    else
      printf("owner_col = %d\n", owner);

    break;
  default:
    printf("Unknown objtype_to\n");
  }

  PLA_Abort("PLA_Reduce_old called, not implemented", __LINE__, __FILE__ );
}
