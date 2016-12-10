/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Copy_local_part_from_matrix_to_msc( PLA_Obj Obj_from, PLA_Obj Msc_to )
{
  PLA_Obj
    from_R  = NULL,       to_R  = NULL,
    from_B1 = NULL,       to_B1 = NULL,
    from_11 = NULL,       to_11 = NULL;

  PLA_Template
    templ;

  int 
    size_left,  size_top, 
    owner_left, owner_top,
    myrow,      mycol;

  PLA_Obj_template( Obj_from, &templ );
  PLA_Temp_comm_col_rank( templ, &myrow );
  PLA_Temp_comm_row_rank( templ, &mycol );

  PLA_Obj_set_to_zero( Msc_to );
    
  PLA_Obj_view_all( Obj_from, &from_R );
  PLA_Obj_view_all( Msc_to,   &to_R );
    
  while ( TRUE ){
    PLA_Obj_split_size( from_R, PLA_SIDE_LEFT, &size_left, &owner_left );
    if ( size_left == 0 ) break;

    PLA_Obj_vert_split_2( from_R, size_left, &from_B1, &from_R );
    PLA_Obj_vert_split_2( to_R,   size_left, &to_B1,   &to_R );

    if ( owner_left == mycol ) {
      while( TRUE ){
	PLA_Obj_split_size( from_B1, PLA_SIDE_TOP, &size_top, &owner_top );
	if ( size_top == 0 ) break;

	PLA_Obj_horz_split_2( from_B1, size_top, &from_11, 
                                                  &from_B1 );

	PLA_Obj_horz_split_2( to_B1,   size_top, &to_11, 
                                                  &to_B1 );

	if ( owner_top == myrow ) PLA_Local_copy( from_11, to_11 );
      }
    }
  }

  PLA_Obj_free( &from_R );    PLA_Obj_free( &to_R );
  PLA_Obj_free( &from_B1 );   PLA_Obj_free( &to_B1 );
  PLA_Obj_free( &from_11 );    PLA_Obj_free( &to_11 );
      
  return PLA_SUCCESS;
}


int PLA_Copy_local_part_from_msc_to_matrix( 
		  PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  PLA_Obj
    from_R  = NULL,       to_R  = NULL,
    from_B1 = NULL,       to_B1 = NULL,
    from_11 = NULL,       to_11 = NULL;

  PLA_Template
    templ;

  int 
    size_left,  size_top, 
    owner_left, owner_top,
    myrow,      mycol;

  PLA_Obj_template( Obj_from, &templ );
  PLA_Temp_comm_col_rank( templ, &myrow );
  PLA_Temp_comm_row_rank( templ, &mycol );

  PLA_Obj_set_to_zero( Obj_to );
    
  PLA_Obj_view_all( Obj_from, &from_R );
  PLA_Obj_view_all( Obj_to,   &to_R );
    
  while ( TRUE ){
    PLA_Obj_split_size( to_R, PLA_SIDE_LEFT, &size_left, &owner_left );
    if ( size_left == 0 ) break;

    PLA_Obj_vert_split_2( from_R, size_left, &from_B1, &from_R );
    PLA_Obj_vert_split_2( to_R,   size_left, &to_B1,   &to_R );

    if ( owner_left == mycol ) {
      while( TRUE ){
	PLA_Obj_split_size( to_B1, PLA_SIDE_TOP, &size_top, &owner_top );
	if ( size_top == 0 ) break;

	PLA_Obj_horz_split_2( from_B1, size_top, &from_11, 
                                                  &from_B1 );

	PLA_Obj_horz_split_2( to_B1,   size_top, &to_11, 
                                                  &to_B1 );

	if ( owner_top == myrow ) PLA_Local_copy( from_11, to_11 );
      }
    }
  }

  PLA_Obj_free( &from_R );    PLA_Obj_free( &to_R );
  PLA_Obj_free( &from_B1 );   PLA_Obj_free( &to_B1 );
  PLA_Obj_free( &from_11 );    PLA_Obj_free( &to_11 );
      
  return PLA_SUCCESS;
}

int PLA_Copy_local_part_from_mv_to_msc( PLA_Obj Obj_from, PLA_Obj Msc_to )
{
  PLA_Obj
    from_B = NULL,       to_B = NULL,
    from_1 = NULL,       to_1 = NULL;

  PLA_Template
    templ;

  int 
    size_top, owner_top, me;

  PLA_Obj_template( Obj_from, &templ );
  PLA_Temp_comm_all_rank( templ, &me );

  PLA_Obj_set_to_zero( Msc_to );
    
  PLA_Obj_view_all( Obj_from, &from_B );
  PLA_Obj_view_all( Msc_to,   &to_B );
    
  while ( TRUE ){
    PLA_Obj_split_size( from_B, PLA_SIDE_TOP, &size_top, &owner_top );
    if ( size_top == 0 ) break;

    PLA_Obj_horz_split_2( from_B, size_top, &from_1, 
                                             &from_B );
    PLA_Obj_horz_split_2( to_B,   size_top, &to_1,   
                                              &to_B );

    if ( owner_top == me ) PLA_Local_copy( from_1, to_1 );
  }

  PLA_Obj_free( &from_B );   PLA_Obj_free( &to_B );
  PLA_Obj_free( &from_1 );    PLA_Obj_free( &to_1 );
      
  return PLA_SUCCESS;
}



