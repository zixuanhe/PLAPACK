/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Obj_split_size_to_next_proc( PLA_Obj A, int side, int *size, int *owner )
{
  int size_cur, owner_cur;
  PLA_Obj Acur = NULL;

  PLA_Obj_view_all( A, &Acur );

  PLA_Obj_split_size( Acur, side, size, owner );
  if ( 0 == size ) return;

  switch (side) {
  case PLA_SIDE_TOP:
    PLA_Obj_horz_split_2( Acur, *size, PLA_DUMMY,
			              &Acur );
    while( TRUE ){
      PLA_Obj_split_size( Acur, side, &size_cur, &owner_cur );
      if ( size_cur == 0 || owner_cur != *owner ) break;
      *size += size_cur;
      PLA_Obj_horz_split_2( Acur, size_cur, PLA_DUMMY,
                                            &Acur );
    }
    break;
  case PLA_SIDE_LEFT:
    PLA_Obj_vert_split_2( Acur, *size, PLA_DUMMY, &Acur );
    while( TRUE ){
      PLA_Obj_split_size( Acur, side, &size_cur, &owner_cur );
      if ( size_cur == 0 || owner_cur != *owner ) break;
      *size += size_cur;
      PLA_Obj_vert_split_2( Acur, size_cur, PLA_DUMMY, &Acur );
    }
    break;
  case PLA_SIDE_BOTTOM:
    PLA_Obj_horz_split_2( Acur, - *size, &Acur,
			               PLA_DUMMY );
    while( TRUE ){
      PLA_Obj_split_size( Acur, side, &size_cur, &owner_cur );
      if ( size_cur == 0 || owner_cur != *owner ) break;
      *size += size_cur;
      PLA_Obj_horz_split_2( Acur, -size_cur, &Acur,
                                            PLA_DUMMY);
    }
    break;
  case PLA_SIDE_RIGHT:
    PLA_Obj_vert_split_2( Acur, - *size, &Acur, PLA_DUMMY );
    while( TRUE ){
      PLA_Obj_split_size( Acur, side, &size_cur, &owner_cur );
      if ( size_cur == 0 || owner_cur != *owner ) break;
      *size += size_cur;
      PLA_Obj_vert_split_2( Acur, -size_cur, &Acur, PLA_DUMMY );
    }
  }
  
  PLA_Obj_free( &Acur );

  return PLA_SUCCESS;
}
