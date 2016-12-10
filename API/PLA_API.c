/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include <stdio.h>

/*
   State of PLAPACK (PLA_STATE_NORMAL or PLA_API_ACTIVE)
*/

static int PLA_state = PLA_STATE_NORMAL;

#define MAX_BLK_SIZE    32

/******************************************************************************/

int PLA_API_begin()

/*--------------------------------------------------------------------------

Purpose : Activates the API active state.

----------------------------------------------------------------------------*/
{
  MPI_Comm 
    base_comm = MPI_COMM_NULL;

  if ( PLA_state != PLA_STATE_NORMAL ) 
    PLA_Abort( "API already active", __LINE__, __FILE__ );

  PLA_state = PLA_STATE_API_ACTIVE;

  /* Initialize put buffers */
  PLA_API_Init_put_buffers( );

  /* Initialize receive buffer */
  PLA_API_Init_recv_buffer( );

  /* Initialize request buffers */
  PLA_API_Init_request_buffers( );

  /* Initialize open objects list */
  PLA_API_Init_open_objects_list( );

  /* Synchronize */
  PLA_Base_comm_1d( &base_comm );
  MPI_Barrier( base_comm );

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_end()

/*--------------------------------------------------------------------------

Purpose : De-activates the API active state.



----------------------------------------------------------------------------*/
{
  MPI_Comm 
    comm = MPI_COMM_NULL;
  
  int 
    nprocs,
    number_open_objects;

  pla2_api_sync_request( FALSE ); 
  pla2_api_sync_put( FALSE ); 

  PLA_Base_comm_1d( &comm );

  if ( PLA_state != PLA_STATE_API_ACTIVE ) 
    PLA_Abort( "API not active", __LINE__, __FILE__ );

  PLA_API_number_open_objects( &number_open_objects );
  if ( number_open_objects != 0 )
    PLA_Abort( "Cannot end API without closing all objects", __LINE__, __FILE__ );

  MPI_Barrier( comm);

  PLA_API_flush_sync_messages( );

  PLA_state = PLA_STATE_NORMAL;

  PLA_API_Free_open_objects_list( );

  PLA_API_Free_request_buffers( );

  PLA_API_Free_recv_buffer( );

  PLA_API_Free_put_buffers( );

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_state(int    *state)

/*--------------------------------------------------------------------------

Purpose : Determines the state of a linear algebra object.

OUT     state    State of linear algebra object

----------------------------------------------------------------------------*/
{
  *state = PLA_state;

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_Obj_API_open(PLA_Obj obj)

/*--------------------------------------------------------------------------

Purpose : Puts object into asynchronous mode
          (makes PLA_Obj_API_ operations on the object allowable).


IN/OUT     obj    Linear algebra object (to be) placed in asynchronous mode.

----------------------------------------------------------------------------*/
{
  PLA_Template
    template = NULL;
  
  MPI_Comm 
    comm = MPI_COMM_NULL;
  
  int 
    mode;

  PLA_Obj_template( obj, &template );
  PLA_Temp_comm_all( template, &comm );

  MPI_Barrier( comm );

  PLA_Obj_API_mode( obj, &mode );

  if ( mode != PLA_MODE_CLOSED )
    PLA_Abort( "object already open", __LINE__, __FILE__ );

  PLA_API_add_obj_to_open_objects_list( obj );

  obj->mode = PLA_MODE_OPEN;

  return PLA_SUCCESS; 
}

/******************************************************************************/

int PLA_Obj_API_close( PLA_Obj obj )

/*--------------------------------------------------------------------------

Purpose : Puts object into synchronous mode (disallows PLA_Obj_API_
operations on the object).


IN/OUT     obj    Linear algebra object (to be) taken out of asynchronous mode.

----------------------------------------------------------------------------*/
{
  int 
    mode,
    i;

  PLA_Obj_API_mode( obj, &mode );

  if ( mode != PLA_MODE_OPEN )
    PLA_Abort( "object not open", __LINE__, __FILE__ );

  /* Synchronize API, reposting requests, since API_end has not been called */
  PLA_Obj_API_sync( obj );

  PLA_API_delete_obj_from_open_objects_list( obj );

  /* Mark object as closed */
  obj->mode = PLA_MODE_CLOSED;

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_Obj_API_mode (PLA_Obj  obj,     int       *mode)

/*--------------------------------------------------------------------------

Purpose : Returns the mode of the linear algebra object (PLA_MODE_CLOSED
or PLA_MODE_OPEN).


IN     obj       Linear algebra object
OUT    mode      Mode of linear algebra object

----------------------------------------------------------------------------*/
{
  *mode = obj->mode; 

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_axpy_vector_to_global(   int            size ,
                                     void         * alpha ,
                                     void         * local_obj ,
                                     int            local_inc ,
                                     PLA_Obj      laobj ,
                                     int            disp )

/*--------------------------------------------------------------------------

Purpose : Add local vector buffer data to global vector object.


IN            size              Length of vector
IN          * alpha             Scalar in axpy operation
IN          * local_obj         Local vector values
IN            local_inc         Local vector increment (stride)
IN/OUT        laobj             Global vector (target linear algebra object)
IN            disp              Global vector displacement

----------------------------------------------------------------------------*/
{
  PLA_Template templ;
  int zero_or_one;

  PLA_Obj_template( laobj, &templ );
  PLA_Temp_zero_or_one( templ, &zero_or_one );
  if ( local_inc == 1 )
    PLA_API_axpy_matrix_to_global( size, 1, alpha, local_obj, size, 
            laobj, disp, zero_or_one );
  else 
    PLA_Abort( "API_axpy_vector_to_global: stride != 1 not implemented", __LINE__, __FILE__ );
    
  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_axpy_global_to_vector(int           size ,
                                  void        * alpha ,
                                  PLA_Obj     laobj ,
                                  int           displ ,
                                  void        * local_obj ,
                                  int           local_inc )

/*--------------------------------------------------------------------------

Purpose : Add piece of global vector object to local linear algebra
object.


IN            size              Length of vector
IN          * alpha             Scalar in axpy operation
IN            laobj             Global vector object
IN            displ             Displacement  (in global object)
IN/OUT      * local_buf         Local linear algebra object (target) buffer
IN            local_inc         Local vector increment (stride)

----------------------------------------------------------------------------*/
{
  PLA_Template templ;
  int zero_or_one;

  PLA_Obj_template( laobj, &templ );
  PLA_Temp_zero_or_one( templ, &zero_or_one );
  if ( local_inc == 1 )
    PLA_API_axpy_global_to_matrix( size, 1, alpha, laobj, displ, zero_or_one,
            local_obj, size );
  else 
    PLA_Abort( "API_axpy_global_to_vector: stride != 1 not implemented", __LINE__, __FILE__ );
    
  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_multi_axpy_vector_to_global(int            nsub,
                                           int          * sizes ,
                                           void         * alpha ,
                                           void         * local_vector ,
                                           int            local_stride ,
                                           PLA_Obj      laobj ,
                                           int          * displs )

/*--------------------------------------------------------------------------

Purpose : Add local sub-vector buffer data chunks to global vector
object.


IN            nsub              Number of subvectors
IN          * size              Length of subvectors
IN          * alpha             Scalar in axpy operation
IN          * local_vector      Local vector values
IN            local_stride      Local vector increment (stride)
IN/OUT        laobj             Global vector (target linear algebra object)
IN          * displs            Global vector displacements

----------------------------------------------------------------------------*/
{
  int
    i,
    typesize,
    local_i;

  MPI_Datatype
    datatype;

  PLA_Obj_datatype( laobj, &datatype );
  MPI_Type_size( datatype, &typesize );

  local_i = 0;
  for ( i=0; i<nsub; i++ ){
    PLA_API_axpy_vector_to_global( sizes[ i ], 
              alpha, (char *) local_vector + local_i*typesize,
              local_stride, laobj, displs[ i ] );
    local_i += sizes[ i ];
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_multi_axpy_global_to_vector(int           nsub ,
                                           int         * sizes ,
                                           void        * alpha ,
                                           PLA_Obj     laobj ,
                                           int         * displs ,
                                           void        * local_vector ,
                                           int           local_stride )

/*--------------------------------------------------------------------------

Purpose : Add pieces of global vector object to local vector.

 IN            nsub              Number of subvectors
IN          * size              Length of subvectors
IN          * alpha             Scalar in axpy operation
IN            laobj             Global vector object
IN          * displs            Displacements  (in global object)
IN/OUT      * local_vector      Local linear algebra object (target) buffer
IN            local_stride      Local vector increment (stride)

----------------------------------------------------------------------------*/
{
  int
    i,
    typesize,
    local_i;

  MPI_Datatype
    datatype;

  PLA_Obj_datatype( laobj, &datatype );
  MPI_Type_size( datatype, &typesize );

  local_i = 0;
  for ( i=0; i<nsub; i++ ){
    PLA_API_axpy_global_to_vector( sizes[ i ], 
              alpha, laobj, displs[ i ], 
              (char *) local_vector + local_i*typesize,
              local_stride );
    local_i += sizes[ i ];
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_axpy_matrix_to_global(   int            m,
                                     int            n,
                                     void         * alpha ,
                                     void         * local_buf, 
                                     int            local_ldim ,
                                     PLA_Obj       A,
                                     int            disp_row,
                                     int            disp_col )

/*--------------------------------------------------------------------------

Purpose : Add local matrix buffer data to global matrix object.


IN            m                 Length of matrix
IN            n                 Width of matrix
IN          * alpha             Scalar in axpy operation
IN          * local_buf         Local matrix values (buffer)
IN            local_ldim        Local matrix leading dimension
IN/OUT        A                 Global matrix (target linear algebra object)
IN            disp_row          Global matrix row displacement
IN            disp_col          Global matrix column displacement

----------------------------------------------------------------------------*/
{
  MPI_Datatype
    datatype;

  MPI_Comm
    comm = MPI_COMM_NULL;

  PLA_Template
    templ = NULL;

  PLA_Obj
    A_R  = NULL, A_B1 = NULL, A_11 = NULL;

  int 
    objtype,
    typesize,
    me, 
    myrow, mycol, nprows, 
    size_left, size_top,
    owner_left, owner_top,
    obj_index,
    local_i, local_j,
    global_i, global_j,
    proc_index;

  PLA_Obj_objtype( A,  &objtype );
  PLA_Obj_datatype( A, &datatype );
  MPI_Type_size( datatype, &typesize );
  PLA_Obj_template( A, &templ );

  PLA_Obj_index_in_open_objects( A, &obj_index );
  
  switch ( objtype ){
  case PLA_MATRIX:
    PLA_Temp_comm_col_rank( templ, &myrow );
    PLA_Temp_comm_col_size( templ, &nprows );
    PLA_Temp_comm_row_rank( templ, &mycol );

    /* A_R views the part of the object to which the data is
       to be added */
    PLA_Obj_view( A, m, n, disp_row, disp_col, &A_R );

    /* local_j and global_j keep track of the column in the
       local data and global object, respectively */
    local_j = 0;
    global_j = disp_col;

    while ( TRUE ){
      /* Determine size and owner of column panel that resides on
         one column of nodes */
      PLA_Obj_split_size( A_R, PLA_SIDE_LEFT, &size_left, &owner_left );
      if ( size_left == 0 ) break;
      size_left = min( size_left, MAX_BLK_SIZE );

      /* Partition A_R = < A_B1, A_R > where A_B1 resides within
         one column of nodes */
      PLA_Obj_vert_split_2( A_R, size_left, &A_B1, &A_R );

      /* local_i and global_i keep track of the row in the
         local data and global object, respectively */
      local_i = 0;
      global_i = disp_row;

      while ( TRUE ){
	/* Determine size and owner of row panel that resides on
           one row of nodes */
	PLA_Obj_split_size( A_B1, PLA_SIDE_TOP, &size_top, &owner_top);
	if ( size_top == 0 ) break;
	size_top = min( size_top, MAX_BLK_SIZE );
      
	/* Partition A_B1 = / A_11 \  where A_11 resides within
                            \ A_B1 /  one node */
	PLA_Obj_horz_split_2( A_B1, size_top, &A_11,
                                              &A_B1 );

	if ( myrow == owner_top && mycol == owner_left ){
	  /* Data needs to be entered locally */
	  PLA_Obj_add_to_local_contents( 
              PLA_NO_TRANS, size_top, size_left,
              ( void * ) ( ( char * )local_buf + 
              ( local_j * local_ldim + local_i ) * typesize ),
	      local_ldim, 1, A_11 );
	}
	else {
	  proc_index = owner_left * nprows + owner_top;
	  PLA_API_put( datatype, size_top, size_left,
              ( void * ) ( ( char * )local_buf + 
              ( local_j * local_ldim + local_i ) * typesize ),
	      local_ldim, proc_index, obj_index, global_i,
	      global_j );
	}
	global_i += size_top;
	local_i += size_top;
      }
      global_j += size_left;
      local_j += size_left;
    }
    
    PLA_Obj_free( &A_R );
    PLA_Obj_free( &A_B1 );
    PLA_Obj_free( &A_11 );
    break;

  case PLA_MVECTOR:
    PLA_Temp_comm_all_rank( templ, &me );

    /* A_R views the part of the object to which the data is
       to be added */
    PLA_Obj_view( A, m, n, disp_row, disp_col, &A_R );

    /* local_j and global_j keep track of the column in the
       local data and global object, respectively */
    local_j = 0;
    global_j = disp_col;

    while ( TRUE ){
      /* Determine size of column panel */
      PLA_Obj_global_width( A_R, &size_left );
      if ( size_left == 0 ) break;
      size_left = min( size_left, MAX_BLK_SIZE );

      /* Partition A_R = < A_B1, A_R > where A_B1 resides within
         one column of nodes */
      PLA_Obj_vert_split_2( A_R, size_left, &A_B1, &A_R );

      /* local_i and global_i keep track of the row in the
         local data and global object, respectively */
      local_i = 0;
      global_i = disp_row;

      while ( TRUE ){
	/* Determine size and owner of row panel that resides on
           one row of nodes */
	PLA_Obj_split_size( A_B1, PLA_SIDE_TOP, &size_top, &owner_top);
	if ( size_top == 0 ) break;
	size_top = min( size_top, MAX_BLK_SIZE );
      
	/* Partition A_B1 = / A_11 \  where A_11 resides within
                            \ A_B1 /  one node */
	PLA_Obj_horz_split_2( A_B1, size_top, &A_11,
                                              &A_B1 );

	if ( me == owner_top ) {
	  /* Data needs to be entered locally */
	  PLA_Obj_add_to_local_contents( 
              PLA_NO_TRANS, size_top, size_left,
              ( void * ) ( ( char * )local_buf + 
              ( local_j * local_ldim + local_i ) * typesize ),
	      local_ldim, 1, A_11 );
	}
	else {
	  proc_index = owner_top;
	  PLA_API_put( datatype, size_top, size_left,
              ( void * ) ( ( char * )local_buf + 
              ( local_j * local_ldim + local_i ) * typesize ),
	      local_ldim, proc_index, obj_index, global_i,
	      global_j );
	}
	global_i += size_top;
	local_i += size_top;
      }
      global_j += size_left;
      local_j += size_left;
    }
    
    PLA_Obj_free( &A_R );
    PLA_Obj_free( &A_B1 );
    PLA_Obj_free( &A_11 );
    break;

  default:
    PLA_Abort( "Illegal target object type", __LINE__, __FILE__ );
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_axpy_global_to_matrix(int           m,
                                  int           n,
                                  void        * alpha ,
                                  PLA_Obj       A,
                                  int           displ_row ,
                                  int           displ_col,
                                  void        * local_buf,
                                  int           local_ldim )

/*--------------------------------------------------------------------------

Purpose : Add piece of global matrix object to local linear algebra
object (matrix).


IN            m                 Length of matrix
IN            n                 Width of matrix
IN          * alpha             Scalar in axpy operation
IN            A                 Global object
IN            displ_row         Row displacement  (in global object)
IN            displ_col         Column displacement  (in global object)
IN/OUT      * local_buf      Local linear algebra object (target) buffer
IN            local_ldim         Local matrix leading dimension

----------------------------------------------------------------------------*/
{
  MPI_Datatype
    datatype;

  MPI_Comm
    comm = MPI_COMM_NULL;

  PLA_Template
    templ = NULL;

  PLA_Obj
    A_R  = NULL, A_B1 = NULL, A_11 = NULL;

  int 
    objtype,
    typesize,
    me, 
    myrow, mycol, nprows, 
    size_left, size_top,
    owner_left, owner_top,
    obj_index,
    local_i, local_j,
    global_i, global_j,
    proc_index;

  PLA_Obj_objtype( A,  &objtype );
  PLA_Obj_datatype( A, &datatype );
  MPI_Type_size( datatype, &typesize );
  PLA_Obj_template( A, &templ );

  PLA_Obj_index_in_open_objects( A, &obj_index );
  
  switch ( objtype ){
  case PLA_MATRIX:
    PLA_Temp_comm_col_rank( templ, &myrow );
    PLA_Temp_comm_col_size( templ, &nprows );
    PLA_Temp_comm_row_rank( templ, &mycol );

    /* A_R views the part of the object to which the data is
       to be added */
    PLA_Obj_view( A, m, n, displ_row, displ_col, &A_R );

    /* local_j and global_j keep track of the column in the
       local data and global object, respectively */
    local_j = 0;
    global_j = displ_col;

    while ( TRUE ){
      /* Determine size and owner of column panel that resides on
         one column of nodes */
      PLA_Obj_split_size( A_R, PLA_SIDE_LEFT, &size_left, &owner_left );
      if ( size_left == 0 ) break;
      size_left = min( size_left, MAX_BLK_SIZE );

      /* Partition A_R = < A_B1, A_R > where A_B1 resides within
         one column of nodes */
      PLA_Obj_vert_split_2( A_R, size_left, &A_B1, &A_R );

      /* local_i and global_i keep track of the row in the
         local data and global object, respectively */
      local_i = 0;
      global_i = displ_row;

      while ( TRUE ){
	/* Determine size and owner of row panel that resides on
           one row of nodes */
	PLA_Obj_split_size( A_B1, PLA_SIDE_TOP, &size_top, &owner_top);
	if ( size_top == 0 ) break;
	size_top = min( size_top, MAX_BLK_SIZE );
      
	/* Partition A_B1 = / A_11 \  where A_11 resides within
                            \ A_B1 /  one node */
	PLA_Obj_horz_split_2( A_B1, size_top, &A_11,
                                              &A_B1 );

	if ( myrow == owner_top && mycol == owner_left ){
	  /* Data needs to be entered locally */
	  PLA_Obj_add_from_local_contents( 
              PLA_NO_TRANS, size_top, size_left,
              ( void * ) ( ( char * )local_buf + 
              ( local_j * local_ldim + local_i ) * typesize ),
	      local_ldim, 1, A_11 );
	}
	else {
	  proc_index = owner_left * nprows + owner_top;
	  PLA_API_request( datatype, size_top, size_left,
              ( void * ) ( ( char * )local_buf + 
              ( local_j * local_ldim + local_i ) * typesize ),
	      local_ldim, proc_index, obj_index, global_i,
	      global_j );
	}
	global_i += size_top;
	local_i += size_top;
      }
      global_j += size_left;
      local_j += size_left;
    }
    
    PLA_Obj_free( &A_R );
    PLA_Obj_free( &A_B1 );
    PLA_Obj_free( &A_11 );
    break;

  case PLA_MVECTOR:
    PLA_Temp_comm_all_rank( templ, &me );

    /* A_R views the part of the object to which the data is
       to be added */
    PLA_Obj_view( A, m, n, displ_row, displ_col, &A_R );

    /* local_j and global_j keep track of the column in the
       local data and global object, respectively */
    local_j = 0;
    global_j = displ_col;

    while ( TRUE ){
      /* Determine size of column panel */
      PLA_Obj_global_width( A_R, &size_left );
      if ( size_left == 0 ) break;
      size_left = min( size_left, MAX_BLK_SIZE );

      /* Partition A_R = < A_B1, A_R > where A_B1 resides within
         one column of nodes */
      PLA_Obj_vert_split_2( A_R, size_left, &A_B1, &A_R );

      /* local_i and global_i keep track of the row in the
         local data and global object, respectively */
      local_i = 0;
      global_i = displ_row;

      while ( TRUE ){
	/* Determine size and owner of row panel that resides on
           one row of nodes */
	PLA_Obj_split_size( A_B1, PLA_SIDE_TOP, &size_top, &owner_top);
	if ( size_top == 0 ) break;
	size_top = min( size_top, MAX_BLK_SIZE );
      
	/* Partition A_B1 = / A_11 \  where A_11 resides within
                            \ A_B1 /  one node */
	PLA_Obj_horz_split_2( A_B1, size_top, &A_11,
                                              &A_B1 );

	if ( me == owner_top ) {
	  /* Data needs to be entered locally */
	  PLA_Obj_add_from_local_contents( 
              PLA_NO_TRANS, size_top, size_left,
              ( void * ) ( ( char * )local_buf + 
              ( local_j * local_ldim + local_i ) * typesize ),
	      local_ldim, 1, A_11 );
	}
	else {
	  proc_index = owner_top;
	  PLA_API_request( datatype, size_top, size_left,
              ( void * ) ( ( char * )local_buf + 
              ( local_j * local_ldim + local_i ) * typesize ),
	      local_ldim, proc_index, obj_index, global_i,
	      global_j );
	}
	global_i += size_top;
	local_i += size_top;
      }
      global_j += size_left;
      local_j += size_left;
    }
    
    PLA_Obj_free( &A_R );
    PLA_Obj_free( &A_B1 );
    PLA_Obj_free( &A_11 );
    break;

  default:
    PLA_Abort( "Illegal target object type", __LINE__, __FILE__ );
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_multi_axpy_matrix_to_global(   int            nsub_row,
                                           int            nsub_col,
                                           int          * size_row ,
                                           int          * size_col ,
                                           void         * alpha ,
                                           void         * local_matrix ,
                                           int            local_ldim ,
                                           PLA_Obj      obj ,
                                           int          * disp_row,
                                           int          * disp_col )

/*--------------------------------------------------------------------------

Purpose : Add local sub-matrix buffer data chunks to global matrix
object.


IN            nsub_row          Number of row subblocks to map
IN            nsub_row          Number of column subblocks to map
IN          * size_row          Row block lengths
IN          * size_col          Column block widths
IN          * alpha             Scalar in axpy operation
IN          * local_matrix      Local matrix values
IN            local_ldim        Local matrix leading dimension
IN/OUT        obj               Global matrix (target linear algebra object)
IN          * disp_row          Global matrix row displacements
IN          * disp_col          Global matrix column displacements

----------------------------------------------------------------------------*/
{
  int
    i, j,
    typesize,
    local_i, local_j;

  MPI_Datatype
    datatype;

  PLA_Obj_datatype( obj, &datatype );
  MPI_Type_size( datatype, &typesize );

  local_i = 0;
  for ( i=0; i<nsub_row; i++ ){
    local_j = 0;
    for ( j=0; j<nsub_col; j++ ){
      PLA_API_axpy_matrix_to_global( size_row[ i ], size_col[ j ],
              alpha, (char *) local_matrix + 
              (local_j*local_ldim + local_i)*typesize,
              local_ldim, obj, disp_row[ i ], disp_col[ j ] );
      local_j += size_col[ j ];
    }
    local_i += size_row[ i ];
  }

  return PLA_SUCCESS;
}

/****************************************************************************/

int PLA_API_multi_axpy_global_to_matrix(   int          nsub_row ,
                                           int           nsub_col ,
                                           int         * size_row ,
                                           int         * size_col ,
                                           void        * alpha ,
                                           PLA_Obj       obj ,
                                           int         * disp_row ,
                                           int         * disp_col ,
                                           void        * local_matrix ,
                                           int           local_ldim )

/*--------------------------------------------------------------------------

Purpose : Add pieces of global matrix object to local matrix.



IN            nsub_row          Number of row subblocks to map
IN            nsub_row          Number of column subblocks to map
IN          * size_row          Row block lengths
IN          * size_col          Column block widths
IN          * alpha             Scalar in axpy operation
IN            obj               Global matrix
IN          * disp_row          Global matrix row displacements
IN          * disp_col          Global matrix column displacements
IN/OUT      * local_matrix      Local matrix values
IN            local_ldim        Local matrix leading dimension

----------------------------------------------------------------------------*/
{
  int
    i, j,
    typesize,
    local_i, local_j;

  MPI_Datatype
    datatype;

  PLA_Obj_datatype( obj, &datatype );
  MPI_Type_size( datatype, &typesize );

  local_i = 0;
  for ( i=0; i<nsub_row; i++ ){
    local_j = 0;
    for ( j=0; j<nsub_col; j++ ){
      PLA_API_axpy_global_to_matrix( size_row[ i ], size_col[ j ],
              alpha, obj, disp_row[ i ], disp_col[ j ],
              (char *) local_matrix + 
              (local_j*local_ldim + local_i)*typesize,
              local_ldim );
      local_j += size_col[ j ];
    }
    local_i += size_row[ i ];
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_Obj_API_sync(PLA_Obj   laobj)

/*--------------------------------------------------------------------------

Purpose : Synchronizes all nodes w.r.t. the given object and completes
all pending "API_axpy" operations on the object.


IN/OUT    laobj      Linear algebra object

----------------------------------------------------------------------------*/
{
  pla2_api_sync_request( TRUE ); 
  pla2_api_sync_put( TRUE );

  return PLA_SUCCESS;
}


