/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Trmv_perform_local_part( int uplo, int trans, int diag,
				  PLA_Obj A, PLA_Obj x_dup, 
				  PLA_Obj y_dup )
{
  int 
    myrow, mycol, 
    owner_top, owner_bottom, owner_left, owner_right,
    size, size_top,  size_bottom,  size_left,  size_right;

  PLA_Obj
    A_cur = NULL, A_11 = NULL, A_21 = NULL, A_12 = NULL,
    x_dup_cur = NULL, x_dup_1 = NULL,
    y_dup_cur = NULL, y_dup_1 = NULL, one = NULL;

  PLA_Template
    templ;

  PLA_Obj_set_to_zero( y_dup );

  PLA_Obj_template ( A, &templ );
  PLA_Temp_comm_col_rank( templ, &myrow );
  PLA_Temp_comm_row_rank( templ, &mycol );

  PLA_Create_constants_conf_to( A, NULL, NULL, &one );

  if ( trans == PLA_NO_TRANSPOSE ){
    if ( uplo == PLA_LOWER_TRIANGULAR ){
      if ( diag == PLA_ZERO_DIAG ){
	PLA_Obj_split_4( A, 1, -1,   PLA_DUMMY, PLA_DUMMY,
                                      &A_cur,     PLA_DUMMY );

	PLA_Obj_vert_split_2( x_dup, -1, &x_dup_cur, PLA_DUMMY );

	PLA_Obj_horz_split_2( y_dup, 1,  PLA_DUMMY,
			              &y_dup_cur );

	diag = PLA_NONUNIT_DIAG;
      }
      else{
	PLA_Obj_view_all( A,     &A_cur );
	PLA_Obj_view_all( x_dup, &x_dup_cur );
	PLA_Obj_view_all( y_dup, &y_dup_cur );
      }
    
      while ( TRUE ){
	PLA_Obj_split_size( A_cur, PLA_SIDE_BOTTOM, &size_bottom, &owner_bottom );
	PLA_Obj_split_size( A_cur, PLA_SIDE_RIGHT,  &size_right,  &owner_right );
	if ( 0 == ( size = min( size_bottom, size_right ) ) ) break;

	PLA_Obj_split_4( A_cur, -size, -size, &A_cur, PLA_DUMMY,
                                               &A_21, &A_11 );

	PLA_Obj_vert_split_2( x_dup_cur, -size, &x_dup_cur, &x_dup_1 );
	PLA_Obj_horz_split_2( y_dup_cur, -size, &y_dup_cur,
                                                 &y_dup_1 );
	if ( myrow == owner_bottom ){
	  if ( mycol == owner_right ) {
	    void
	      *buff;
	    int 
	      i_one;

	    PLA_Obj_local_buffer( y_dup_1, &buff );
	    PLA_Obj_get_local_contents( x_dup_1, PLA_NO_TRANSPOSE,
					 &i_one, &size, buff, 1, 1 ); 
	    
	    PLA_Local_trmv( uplo, trans, diag, A_11, y_dup_1 );
	  }
	  PLA_Local_gemv( PLA_NO_TRANSPOSE, one, A_21, x_dup_cur, one, y_dup_1 );
	}
      }
    }
    else { /* trans = PLA_UPPER_TRIANGULAR */
      if ( diag == PLA_ZERO_DIAG ){
	PLA_Obj_split_4( A, -1, 1,   PLA_DUMMY, &A_cur,
                                     PLA_DUMMY, PLA_DUMMY );

	PLA_Obj_vert_split_2( x_dup, 1, PLA_DUMMY, &x_dup_cur );

	PLA_Obj_horz_split_2( y_dup, -1,  &y_dup_cur,
			                  PLA_DUMMY );

	diag = PLA_NONUNIT_DIAG;
      }
      else{
	PLA_Obj_view_all( A,     &A_cur );
	PLA_Obj_view_all( x_dup, &x_dup_cur );
	PLA_Obj_view_all( y_dup, &y_dup_cur );
      }
    
      while ( TRUE ){
	PLA_Obj_split_size( A_cur, PLA_SIDE_TOP,   &size_top,  &owner_top );
	PLA_Obj_split_size( A_cur, PLA_SIDE_LEFT,  &size_left, &owner_left );
	if ( 0 == ( size = min( size_top, size_left ) ) ) break;

	PLA_Obj_split_4( A_cur, size, size,  &A_11,    &A_12,
                                             PLA_DUMMY, &A_cur );

	PLA_Obj_vert_split_2( x_dup_cur, size, &x_dup_1, &x_dup_cur );
	PLA_Obj_horz_split_2( y_dup_cur, size, &y_dup_1,
                                               &y_dup_cur );
	if ( myrow == owner_top ){
	  if ( mycol == owner_left ) {
	    void
	      *buff;
	    int 
	      i_one, ldim;

	    PLA_Obj_local_buffer( y_dup_1, &buff );
	    PLA_Obj_local_ldim  ( y_dup_1, &ldim );

	    PLA_Obj_get_local_contents( x_dup_1, PLA_TRANSPOSE,
					 &i_one, &size, buff, ldim, 1 ); 
	    
	    PLA_Local_trmv( uplo, trans, diag, A_11, y_dup_1 );
	  }
	  PLA_Local_gemv( PLA_NO_TRANSPOSE, one, A_12, x_dup_cur, one, y_dup_1 );
	}
      }
    }
  }
  else if ( trans == PLA_TRANSPOSE || trans == PLA_CONJ_TRANS ){
    if ( uplo == PLA_LOWER_TRIANGULAR ){
      if ( diag == PLA_ZERO_DIAG ){
	PLA_Obj_split_4( A, 1, -1,   PLA_DUMMY, PLA_DUMMY,
                                      &A_cur,     PLA_DUMMY );

	PLA_Obj_vert_split_2( y_dup, -1,  &y_dup_cur, PLA_DUMMY );

	PLA_Obj_horz_split_2( x_dup, 1,   PLA_DUMMY,
			                   &x_dup_cur ); 
	diag = PLA_NONUNIT_DIAG;
      }
      else{
	PLA_Obj_view_all( A,     &A_cur );
	PLA_Obj_view_all( x_dup, &x_dup_cur );
	PLA_Obj_view_all( y_dup, &y_dup_cur );
      }
    
      while ( TRUE ){
	PLA_Obj_split_size( A_cur, PLA_SIDE_TOP,  &size_top, &owner_top );
	PLA_Obj_split_size( A_cur, PLA_SIDE_LEFT, &size_left,&owner_left );
	if ( 0 == ( size = min( size_top, size_left ) ) ) break;

	PLA_Obj_split_4( A_cur, size, size, &A_11, PLA_DUMMY,
                                            &A_21, &A_cur );

	PLA_Obj_vert_split_2( y_dup_cur, size, &y_dup_1, &y_dup_cur );
	PLA_Obj_horz_split_2( x_dup_cur, size, &x_dup_1,
			                       &x_dup_cur );

	if ( mycol == owner_left) {
	  if ( myrow == owner_top ){
	    void
	      *buff;
	    int 
	      i_one, ldim;
	    PLA_Obj_local_buffer( y_dup_1, &buff );
	    PLA_Obj_local_ldim  ( y_dup_1, &ldim );
	    
	    PLA_Obj_get_local_contents( x_dup_1, PLA_TRANSPOSE,
					 &i_one, &size, buff, ldim, 1 ); 

	    PLA_Local_trmv( uplo, trans, diag, A_11, y_dup_1 );
	  }
	  PLA_Local_gemv( trans, one, A_21, x_dup_cur, one, y_dup_1 );
	}
      }
    }
    else { /* uplo = PLA_UPPER_TRIANGULAR */
      if ( diag == PLA_ZERO_DIAG ){
	PLA_Obj_split_4( A, -1, 1,   PLA_DUMMY, &A_cur,
                                     PLA_DUMMY, PLA_DUMMY );

	PLA_Obj_vert_split_2( x_dup, 1, PLA_DUMMY, &x_dup_cur );

	PLA_Obj_horz_split_2( y_dup, -1,  &y_dup_cur,
			                  PLA_DUMMY );

	diag = PLA_NONUNIT_DIAG;
      }
      else{
	PLA_Obj_view_all( A,     &A_cur );
	PLA_Obj_view_all( x_dup, &x_dup_cur );
	PLA_Obj_view_all( y_dup, &y_dup_cur );
      }
    
      while ( TRUE ){
	PLA_Obj_split_size( A_cur, PLA_SIDE_BOTTOM, &size_bottom,  &owner_bottom );
	PLA_Obj_split_size( A_cur, PLA_SIDE_RIGHT,  &size_right, &owner_right );
	if ( 0 == ( size = min( size_bottom, size_right ) ) ) break;

	PLA_Obj_split_4( A_cur, -size, -size,  &A_cur,    &A_12,
                                               PLA_DUMMY, &A_11 );

	PLA_Obj_vert_split_2( y_dup_cur, -size, &y_dup_cur, &y_dup_1 );
	PLA_Obj_horz_split_2( x_dup_cur, -size, &x_dup_cur,
                                                &x_dup_1 );
	if ( mycol == owner_right ){
	  if ( myrow == owner_bottom ) {
	    void
	      *buff;
	    int 
	      i_one, ldim;

	    PLA_Obj_local_buffer( y_dup_1, &buff );
	    PLA_Obj_local_ldim  ( y_dup_1, &ldim );

	    PLA_Obj_get_local_contents( x_dup_1, PLA_TRANSPOSE,
					 &i_one, &size, buff, ldim, 1 ); 
	    
	    PLA_Local_trmv( uplo, trans, diag, A_11, y_dup_1 );
	  }
	  PLA_Local_gemv( trans, one, A_12, x_dup_cur, one, y_dup_1 );
	}
      }
    }
  }
  else 
    PLA_Abort( "case not yet implemented", __LINE__, __FILE__ );

  PLA_Obj_free( &A_cur );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_21 );
  PLA_Obj_free( &A_12 );
  PLA_Obj_free( &x_dup_cur );
  PLA_Obj_free( &x_dup_1 );
  PLA_Obj_free( &y_dup_cur );
  PLA_Obj_free( &y_dup_1 );
  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}
	    




