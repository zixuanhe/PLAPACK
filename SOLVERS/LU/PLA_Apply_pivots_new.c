/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

PLA_Apply_pivots_to_rows( PLA_Obj A, PLA_Obj ipiv )
{
  int 
    objtype_A, objtype_ipiv,
    owner_row, owner_col,
    owner_top, owner_pivot, me,
    size, width, i_one, nb_alg = 128, index, *ipiv_1_p,
    typesize, dummy;

  PLA_Obj
    A_B = NULL, A_B_rest = NULL, 
    a_top = NULL, a_pivot = NULL,
    ipiv_B = NULL, ipiv_1 = NULL, ipiv_1_dup = NULL;

  PLA_Template
    templ;

  MPI_Comm
    comm = MPI_COMM_NULL;

  MPI_Datatype
    datatype;

  void *
    *buffer_from, *buffer_to;

  MPI_Status
    status;

  MPI_Request
    request;

  PLA_Obj_local_width( A, &width );

  PLA_Obj_objtype( A, &objtype_A );
  PLA_Obj_objtype( ipiv, &objtype_ipiv );

  if ( objtype_ipiv == PLA_MVECTOR ){
    PLA_Obj_view_all( A, &A_B );
    PLA_Obj_view_all( ipiv, &ipiv_B );
      
    while ( TRUE ){
      PLA_Obj_global_length( ipiv_B, &size );
      if ( 0 == ( size = min( size, nb_alg ) ) ) break;
      
      PLA_Obj_horz_split_2( ipiv_B, size, &ipiv_1,
			                   &ipiv_B );

      PLA_Mscalar_create_conf_to( ipiv_1, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &ipiv_1_dup );

      PLA_Copy( ipiv_1, ipiv_1_dup );

      PLA_Apply_pivots_to_rows( A_B, ipiv_1_dup );

      PLA_Obj_horz_split_2( A_B, size, PLA_DUMMY,
			                &A_B );
    }
  }
  else if ( objtype_ipiv == PLA_MSCALAR ){
    PLA_Obj_owner_row( ipiv, &owner_row );
    PLA_Obj_owner_col( ipiv, &owner_col );

    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( ipiv, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &ipiv_B );

      PLA_Copy( ipiv, ipiv_B );

      PLA_Apply_pivots_to_rows( A_B, ipiv_B );
    }
    else{
      PLA_Obj_view_all( A, &A_B );
      PLA_Obj_view_all( ipiv, &ipiv_B );
    
      PLA_Obj_template ( A, &templ );
      if ( objtype_A == PLA_MVECTOR ){
	PLA_Temp_comm_all( templ, &comm );
	PLA_Temp_comm_all_rank( templ, &me );
      }
      else /* if ( objtype_A == PLA_MATRIX ) */{
	PLA_Temp_comm_col( templ, &comm );
	PLA_Temp_comm_col_rank( templ, &me );
      }
	
      PLA_Obj_datatype( A, &datatype );
      MPI_Type_size( datatype, &typesize );

      buffer_from = PLA_malloc( typesize * width );
      buffer_to   = PLA_malloc( typesize * width );

      while ( TRUE ){
	PLA_Obj_global_length( ipiv_B, &size );
	if ( size == 0 ) break;
	
	PLA_Obj_horz_split_2( ipiv_B, 1, &ipiv_1,
			                  &ipiv_B );

	PLA_Obj_local_buffer( ipiv_1, (void **) &ipiv_1_p );

	index = *ipiv_1_p;
	if ( index != 0 ){
	  PLA_Obj_horz_split_2( A_B, 1, &a_top,
				         PLA_DUMMY );
	  
	  PLA_Obj_horz_split_2( A_B, index, PLA_DUMMY,
				             &A_B_rest );
	  
	  PLA_Obj_horz_split_2( A_B_rest, 1, &a_pivot,
				              PLA_DUMMY );

	  PLA_Obj_split_size( a_top,   PLA_SIDE_TOP, &dummy, &owner_top );
	  PLA_Obj_split_size( a_pivot, PLA_SIDE_TOP, &dummy, &owner_pivot );

	  if ( owner_top == owner_pivot ){
	    if ( me == owner_top ){
	      PLA_Obj_get_local_contents( a_top, PLA_NO_TRANSPOSE,
					   &i_one, &width, buffer_from, 1, 1 );
	      PLA_Local_copy( a_pivot, a_top );
	      PLA_Obj_set_local_contents( PLA_NO_TRANSPOSE,
					   1, width, buffer_from, 1, 1,
					   a_pivot );
	    }
	  }
	  else {
	    if ( owner_top == me ){
	      PLA_MPI_Irecv( buffer_to, width, datatype, owner_pivot, 9999, comm, 
		        &request );
	      PLA_Obj_get_local_contents( a_top, PLA_NO_TRANSPOSE,
					   &i_one, &width, buffer_from, 1, 1 );
	      PLA_MPI_Send( buffer_from, width, datatype, owner_pivot, 9999, comm );
	      PLA_MPI_Wait( &request, &status );
	      PLA_Obj_set_local_contents( PLA_NO_TRANSPOSE, 1, width, 
					   buffer_to, 1, 1, a_top );
	    }

	    if ( owner_pivot == me ){
	      PLA_MPI_Irecv( buffer_to, width, datatype, owner_top, 9999, comm,
		        &request );
	      PLA_Obj_get_local_contents( a_pivot, PLA_NO_TRANSPOSE,
					   &i_one, &width, buffer_from, 1, 1 );
	      PLA_MPI_Send( buffer_from, width, datatype, owner_top, 9999, comm );
	      PLA_MPI_Wait( &request, &status );
	      PLA_Obj_set_local_contents( PLA_NO_TRANSPOSE, 1, width, 
					   buffer_to, 1, 1, a_pivot );
	    }
	  }
	}
	
	PLA_Obj_horz_split_2( A_B, 1, PLA_DUMMY,
			               &A_B );
      }
      PLA_free( buffer_from );
      PLA_free( buffer_to );
    }
  }

  PLA_Obj_free( &A_B );
  PLA_Obj_free( &A_B_rest );
  PLA_Obj_free( &ipiv_B );
  PLA_Obj_free( &ipiv_1 );
  PLA_Obj_free( &ipiv_1_dup );
  PLA_Obj_free( &a_top );
  PLA_Obj_free( &a_pivot );

  return PLA_SUCCESS;
}






/************************************************************
  Routine to pivot rows of a multivector.

  The top row of the multivector is exchanged with row
  number 'index', and a copy of the new top row is stored
  in row_copy.  If index is 0, no exchange is necessary,
  and the only action is to copy the top row into the
  multiscalar row_copy.

  INPUTS :
  mvector........the object (mxn) containing the rows to be exchanged.
  index..........a 1x1 duplicated mscalar of type integer
                 whose value gives the index of the new top row.

  OUTPUTS :
  mvector........contains the same data as before, with the top
                 and index_th rows switched.
  pivot_row_copy.......a duplicated mscalar conformal to one row of the 
                 mvector (i.e. 1xn), on exit it contains a copy 
		 of the new top row, unless it is passed in as
		 PLA_NULL, in which case it is not used.
************************************************************/
int Pivot_rows_mv ( PLA_Obj source, PLA_Obj index, PLA_Obj pivot_row_out)
{

  int 
    return_value = PLA_SUCCESS,
    zero_or_one,
    index_value = 0,
    *index_value_p;
  
  PLA_Obj
    first_row = NULL,
    pivot_row = NULL,
    pivot_row_copy = NULL,
    first_row_copy = NULL,
    mv_remain = NULL;
  PLA_Template templ = NULL;
  MPI_Comm 
    comm_all  = MPI_COMM_NULL;
  MPI_Status 
    status; 
  MPI_Datatype
    datatype;
  int 
    me, owner_first, owner_pivot, width, dummy, n;
  double 
    *buffer_first, *buffer_pivot;

  PLA_Obj_datatype( source, &datatype );

  PLA_Obj_template( source, &templ );
  PLA_Temp_comm_all( templ, &comm_all );
  PLA_Temp_comm_all_rank( templ, &me );
  PLA_Temp_zero_or_one( templ, &zero_or_one );

  PLA_Obj_local_buffer( index, ( void ** ) &index_value_p );
  index_value = *index_value_p;

  PLA_Obj_horz_split_2 (source, 1,       &first_row, 
                                          PLA_DUMMY );

  PLA_Obj_horz_split_2( source, index_value - zero_or_one,  PLA_DUMMY,
                                                            &mv_remain );
  PLA_Obj_horz_split_2( mv_remain, 1,            &pivot_row,
                                                 PLA_DUMMY );
  PLA_Mscalar_create_conf_to (first_row,
	                      PLA_ALL_ROWS, PLA_ALL_COLS,
                              &first_row_copy);

  PLA_Mscalar_create_conf_to (pivot_row,
                              PLA_ALL_ROWS, PLA_ALL_COLS,
			      &pivot_row_copy);

  PLA_Obj_local_buffer( first_row_copy, (void **) &buffer_first );
  PLA_Obj_local_buffer( pivot_row_copy, (void **) &buffer_pivot );

  PLA_Obj_local_width( first_row_copy, &width );

  PLA_Obj_split_size( first_row, PLA_SIDE_TOP,  &dummy, &owner_first );
  PLA_Obj_split_size( pivot_row, PLA_SIDE_TOP,  &dummy, &owner_pivot );

  if ( me == owner_pivot ) {
    PLA_Local_copy ( pivot_row, pivot_row_copy);
  }

  PLA_MPI_Bcast( buffer_pivot, width, datatype, owner_pivot, comm_all );

  PLA_Local_copy ( pivot_row_copy, pivot_row_out );

  if ( owner_pivot == owner_first ) {
    if ( me == owner_first ) {
      PLA_Local_copy ( first_row, pivot_row );
    }
  }
  else {
    if ( me == owner_first ) {
      PLA_Local_copy( first_row, first_row_copy );
      PLA_MPI_Send( buffer_first, width, datatype, owner_pivot, 9999, comm_all );
    }

    if ( me == owner_pivot ) {
      PLA_MPI_Recv( buffer_first, width, datatype, owner_first, 9999,
		   comm_all, &status );
      PLA_Local_copy( first_row_copy, pivot_row );
    }
  }

  if ( me == owner_first ){
    PLA_Local_copy( pivot_row_copy, first_row );
  }

  PLA_Obj_free( &first_row_copy);
  PLA_Obj_free( &first_row);
  PLA_Obj_free( &pivot_row);
  PLA_Obj_free( &pivot_row_copy);
  PLA_Obj_free( &mv_remain );

  return return_value;
}



PLA_Apply_pivots_to_columns_in_reverse( PLA_Obj A, PLA_Obj ipiv )
{
  int 
    objtype_A, objtype_ipiv,
    owner_row, owner_col,
    owner_cur, owner_pivot, me,
    size, length, global_width, i_one, nb_alg = 128, index, *ipiv_cur_p,
    typesize, dummy;

  PLA_Obj
    A_L = NULL, A_R = NULL, A_R_rest = NULL, 
    a_cur = NULL, a_pivot = NULL,
    ipiv_T = NULL, ipiv_1 = NULL, ipiv_1_dup = NULL, ipiv_cur = NULL;

  PLA_Template
    templ;

  MPI_Comm
    comm = MPI_COMM_NULL;

  MPI_Datatype
    datatype;

  void *
    *buffer_from, *buffer_to;

  MPI_Status
    status;

  MPI_Request
    request;

  PLA_Obj_local_length( A, &length );
  PLA_Obj_global_width( A, &global_width );
  
  PLA_Obj_objtype( A, &objtype_A );
  PLA_Obj_objtype( ipiv, &objtype_ipiv );

  if ( objtype_ipiv == PLA_MVECTOR ){
    PLA_Obj_vert_split_2( A, global_width, &A_L, &A_R );
    PLA_Obj_horz_split_2( ipiv, global_width, &ipiv_T,
			                       PLA_DUMMY );
      
    PLA_Obj_template ( A, &templ );
    PLA_Temp_comm_row( templ, &comm );
    PLA_Temp_comm_row_rank( templ, &me );

    PLA_Obj_datatype( A, &datatype );
    MPI_Type_size( datatype, &typesize );

    buffer_from = PLA_malloc( typesize * length );
    buffer_to   = PLA_malloc( typesize * length );

    while ( TRUE ){
      PLA_Obj_global_length( ipiv_T, &size );
      if ( 0 == ( size = min( size, nb_alg ) ) ) break;
      
      PLA_Obj_horz_split_2( ipiv_T, -size, &ipiv_T,
			                    &ipiv_1 );

      PLA_Mscalar_create_conf_to( ipiv_1, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &ipiv_1_dup );

      PLA_Copy( ipiv_1, ipiv_1_dup );

      while ( TRUE ){
	PLA_Obj_global_length( ipiv_1_dup, &size );
	if ( size == 0 ) break;

	PLA_Obj_horz_split_2( ipiv_1_dup, -1, &ipiv_1_dup,
                                               &ipiv_cur );

	PLA_Obj_vert_split_2( A_L, -1, &A_L, &a_cur );

	PLA_Obj_local_buffer( ipiv_cur, (void **) &ipiv_cur_p );

	index = *ipiv_cur_p;
	if ( index != 0 ){
	  PLA_Obj_vert_split_2( A_R, index-1, PLA_DUMMY, &A_R_rest );

	  PLA_Obj_vert_split_2( A_R_rest, 1, &a_pivot, PLA_DUMMY );

	  PLA_Obj_split_size( a_cur,   PLA_SIDE_LEFT, &dummy, &owner_cur );
	  PLA_Obj_split_size( a_pivot, PLA_SIDE_LEFT, &dummy, &owner_pivot );

	  if ( owner_cur == owner_pivot ){
	    if ( me == owner_cur ){
	      PLA_Obj_get_local_contents( a_cur, PLA_NO_TRANSPOSE,
					   &length, &i_one, buffer_from, length, 1 );
	      PLA_Local_copy( a_pivot, a_cur );
	      PLA_Obj_set_local_contents( PLA_NO_TRANSPOSE,
					   length, 1, buffer_from, length, 1,
					   a_pivot );
	    }
	  }
	  else {
	    if ( owner_cur == me ){
	      PLA_MPI_Irecv( buffer_to, length, datatype, owner_pivot, 9999, comm, 
		        &request );
	      PLA_Obj_get_local_contents( a_cur, PLA_NO_TRANSPOSE,
					   &length, &i_one, buffer_from, length, 1 );
	      PLA_MPI_Send( buffer_from, length, datatype, owner_pivot, 9999, comm );
	      PLA_MPI_Wait( &request, &status );
	      PLA_Obj_set_local_contents( PLA_NO_TRANSPOSE, length, 1,
					   buffer_to, length, 1, a_cur );
	    }

	    if ( owner_pivot == me ){
	      PLA_MPI_Irecv( buffer_to, length, datatype, owner_cur, 9999, comm,
		        &request );
	      PLA_Obj_get_local_contents( a_pivot, PLA_NO_TRANSPOSE,
					   &length, &i_one, buffer_from, length, 1 );
	      PLA_MPI_Send( buffer_from, length, datatype, owner_cur, 9999, comm );
	      PLA_MPI_Wait( &request, &status );
	      PLA_Obj_set_local_contents( PLA_NO_TRANSPOSE, length, 1,
					   buffer_to, length, 1, a_pivot );
	    }
	  }
	}
	
	PLA_Obj_view_shift( A_R,    0,
                                 -1,      0,
                                     0 );
      }
    }

    PLA_free( buffer_from );
    PLA_free( buffer_to );
  }
  else
    PLA_Abort("not yet implemented", __LINE__, __FILE__ );

  PLA_Obj_free( &A_R );
  PLA_Obj_free( &A_L );
  PLA_Obj_free( &A_R_rest );
  PLA_Obj_free( &a_cur );
  PLA_Obj_free( &a_pivot );
  PLA_Obj_free( &ipiv_1 );
  PLA_Obj_free( &ipiv_1_dup );
  PLA_Obj_free( &ipiv_cur );

  return PLA_SUCCESS;
}
