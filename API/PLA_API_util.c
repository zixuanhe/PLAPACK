/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include <stdio.h>

#define SIZE_OF_UNIT ( sizeof( double ) )

#define CHECK_FOR_STALL 0
#define STALL_INTERVAL 1.0

/* 
   This implementation of the Application Interface keeps a list
   of open objects, the open_objects_list.
*/

#define MAX_OPEN_OBJS 25

PLA_Obj *open_objects_list;
static int number_open_objects = 0;

struct Send_buffer_struct{
  int
    cur_length;
  char
    *buffer;
  MPI_Request 
    request;
};

/* 
   This implementation of the Application Interface creates
   a "put_buffer" for each processor, in which data to be sent
   from this processor to the other processor is put.
   A list of this information is created, in put_buffer_list.
*/

#define MAX_PUT_BUFFER_SIZE 20000

#define MAX_REQUEST_BUFFER_SIZE 100

static struct Send_buffer_struct 
   *put_buffer_list,
   *request_buffer_list;

/*
   A "put_buffer" from another processor is received into
   the recv_put_buffer;  A current pointer to where in the
   recv buffer a read can occur (much like a file being read)
   is given in recv_put_buffer_curp
*/

/* static char recv_put_buffer[ MAX_PUT_BUFFER_SIZE ]; */
static char *recv_put_buffer;
static int recv_put_buffer_length, recv_put_buffer_done;

static char recv_request_buffer[ MAX_REQUEST_BUFFER_SIZE ]; 
static int recv_request_buffer_length, recv_request_buffer_done;

/*
   Message tags used by the API
*/

#define PUT_READY_TAG   999
#define PUT_TAG   998

#define REQUEST_READY_TAG   997
#define REQUEST_TAG   996

/*
   Types of messages being sent in put_buffer
*/

#define PLA_API_TYPE_SYNC 9999
#define PLA_API_TYPE_PUT  9998
#define PLA_API_TYPE_GET  9997

/*
   Dummy variable into which to receive acknowledgement message
*/

static int dummy;
 
/******************************************************************************/

int PLA_API_Init_open_objects_list( )

/*--------------------------------------------------------------------------

Purpose : Initialize open objects list

----------------------------------------------------------------------------*/
{
  int i;

  open_objects_list = ( PLA_Obj * ) 
    PLA_malloc( (size_t) sizeof( PLA_Obj ) * MAX_OPEN_OBJS );

  if ( open_objects_list == NULL )
    PLA_Abort( "malloc failed", __LINE__, __FILE__ );

  for ( i=0; i<MAX_OPEN_OBJS; i++ ) 
    open_objects_list[ i ] = NULL;
 
  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_Init_put_buffers( )

/*--------------------------------------------------------------------------

Purpose : Initialize put buffers

----------------------------------------------------------------------------*/
{
  int 
    i,
    me, nprocs;

  MPI_Comm 
    base_comm = MPI_COMM_NULL;

  PLA_Base_comm_1d( &base_comm );

  MPI_Comm_size( base_comm, &nprocs );
  MPI_Comm_rank( base_comm, &me );

  put_buffer_list = ( struct Send_buffer_struct * )
    PLA_malloc( (size_t) sizeof( struct Send_buffer_struct ) * nprocs );

  for ( i=0; i<nprocs; i++ ){
    put_buffer_list[ i ].cur_length = 0;

    put_buffer_list[ i ].buffer = 
	( char * ) PLA_malloc( (size_t) MAX_PUT_BUFFER_SIZE );
    if ( put_buffer_list[ i ].buffer == NULL )
      PLA_Abort( "malloc failed", __LINE__, __FILE__ );

    if ( me != i ) {
      MPI_Irecv( &dummy, 1, MPI_INT, i, PUT_READY_TAG, base_comm, 
		&put_buffer_list[ i ].request );

      /* Send signal that this node is ready to receive asynchronouse
	 buffer */

      MPI_Send( &dummy, 1, MPI_INT, i, PUT_READY_TAG, base_comm );
    }
  }

  return PLA_SUCCESS;
}

/************************************************************************/

int PLA_API_Init_recv_buffer( )

/*--------------------------------------------------------------------------

Purpose : Initialize receive buffer

----------------------------------------------------------------------------*/
{
  recv_put_buffer = ( char * ) PLA_malloc( (size_t) MAX_PUT_BUFFER_SIZE );
  if ( recv_put_buffer == NULL )
    PLA_Abort( "malloc failed creating recv_put_buffer", __LINE__, __FILE__ );

  return PLA_SUCCESS;
}

/************************************************************************/

int PLA_API_Init_request_buffers( )

/*--------------------------------------------------------------------------

Purpose : Initialize request buffers

----------------------------------------------------------------------------*/
{
  int 
    i,
    me, nprocs;

  MPI_Comm 
    base_comm = MPI_COMM_NULL;

  PLA_Base_comm_1d( &base_comm );

  MPI_Comm_size( base_comm, &nprocs );
  MPI_Comm_rank( base_comm, &me );

  request_buffer_list = ( struct Send_buffer_struct * )
    PLA_malloc( (size_t) sizeof( struct Send_buffer_struct ) * nprocs );

  for ( i=0; i<nprocs; i++ ){
    request_buffer_list[ i ].cur_length = 0;

    request_buffer_list[ i ].buffer = 
	( char * ) PLA_malloc( (size_t) MAX_REQUEST_BUFFER_SIZE );
    if ( request_buffer_list[ i ].buffer == NULL )
      PLA_Abort( "malloc failed", __LINE__, __FILE__ );

    if ( me != i ) {
      MPI_Irecv( &dummy, 1, MPI_INT, i, REQUEST_READY_TAG, base_comm, 
		&request_buffer_list[ i ].request );

      /* Send signal that this node is ready to receive asynchronouse
	 buffer */

      MPI_Send( &dummy, 1, MPI_INT, i, REQUEST_READY_TAG, base_comm );
    }
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_Free_open_objects_list( )

/*--------------------------------------------------------------------------

Purpose : Free open objects list

----------------------------------------------------------------------------*/
{
  if ( open_objects_list != NULL )
    PLA_free( open_objects_list );

 
  return PLA_SUCCESS;
}

/******************************************************************************/

int   PLA_API_Free_put_buffers( )

/*--------------------------------------------------------------------------

Purpose : Free put buffers

----------------------------------------------------------------------------*/
{
  int 
    i,
    nprocs;

  MPI_Comm 
    comm = MPI_COMM_NULL;

  PLA_Base_comm_1d( &comm );

  MPI_Comm_size( comm, &nprocs );

  if ( put_buffer_list != NULL ) {
    for ( i=0; i<nprocs; i++ )
      if ( put_buffer_list[ i ].buffer != NULL ) 
	PLA_free( put_buffer_list[ i ].buffer );

    PLA_free( put_buffer_list );
  }

  return PLA_SUCCESS;
}


int   PLA_API_Free_recv_buffer( )

/*--------------------------------------------------------------------------

Purpose : Free receive buffer

----------------------------------------------------------------------------*/
{
  if ( recv_put_buffer != NULL ) 
    PLA_free( recv_put_buffer );

  return PLA_SUCCESS;
}

/******************************************************************************/

int   PLA_API_Free_request_buffers( )

/*--------------------------------------------------------------------------

Purpose : Free request buffers

----------------------------------------------------------------------------*/
{
  int 
    i,
    nprocs;

  MPI_Comm 
    comm = MPI_COMM_NULL;

  PLA_Base_comm_1d( &comm );

  MPI_Comm_size( comm, &nprocs );

  if ( request_buffer_list != NULL ) {
    for ( i=0; i<nprocs; i++ )
      if ( request_buffer_list[ i ].buffer != NULL ) 
	PLA_free( request_buffer_list[ i ].buffer );

    PLA_free( request_buffer_list );
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_request  ( MPI_Datatype datatype, 
                       int m, int n, void *local_buf, int lda,
                       int dest, int obj_index, int align_row, int align_col )

/*--------------------------------------------------------------------------

Purpose : enter a request for data on the appropriate request buffer 
          to be sent to processor "dest"

----------------------------------------------------------------------------*/
{
  int 
    me,
    i, j,
    typesize,
    ready,
    temp;

  MPI_Status status;

  MPI_Comm
    base_comm = MPI_COMM_NULL;

  PLA_Base_comm_1d( &base_comm );

  MPI_Comm_rank( base_comm, &me );

  /* Create room to send current message */
  MPI_Type_size( datatype, &typesize );
  PLA_API_request_buffer_make_room( 8 * sizeof( int ) + sizeof( void * ), dest );
			   
  /* Enter information in the request buffer */
  temp = PLA_API_TYPE_GET;
  PLA_API_add_data_to_request_buffer( sizeof( int ), ( char * ) &temp, 
				      dest );
  PLA_API_add_data_to_request_buffer   ( sizeof( int ), ( char * ) &m, 
				       dest );
  PLA_API_add_data_to_request_buffer( sizeof( int ), ( char * ) &n, 
				      dest );
  PLA_API_add_data_to_request_buffer( sizeof( int ), ( char * ) &obj_index, 
				      dest );
  PLA_API_add_data_to_request_buffer( sizeof( int ), ( char * ) &align_row, 
				      dest );
  PLA_API_add_data_to_request_buffer( sizeof( int ), ( char * ) &align_col, 
				      dest );
  PLA_API_add_data_to_request_buffer( sizeof( void * ), ( char * )&local_buf,
				      dest );
  PLA_API_add_data_to_request_buffer( sizeof( int ), ( char * ) &lda,
				      dest );
  PLA_API_add_data_to_request_buffer( sizeof( int ), ( char * ) &me,
				      dest );

  /* Check if anything has come in for this node */
  PLA_API_Poll( );

  return PLA_SUCCESS;
}


int PLA_API_put  ( MPI_Datatype datatype, 
                    int m, int n, void *local_buf, int lda,
                    int dest, int obj_index, int align_row, int align_col )

/*--------------------------------------------------------------------------

Purpose : enter data in the appropriate put buffer to be sent to 
          processor "dest"

----------------------------------------------------------------------------*/
{
  int 
    i, j,
    typesize,
    ready,
    temp;

  char 
    *tempp_to,
    *tempp_from;
    
  MPI_Status status;

  /* Create room to send current message */
  MPI_Type_size( datatype, &typesize );
  PLA_API_put_buffer_make_room( 6 * sizeof( int ) + m * n * typesize, dest );
			   
  /* Enter put destination information in the buffer */
  temp = PLA_API_TYPE_PUT;
  PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &temp,      dest );
  PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &m,         dest );
  PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &n,         dest );
  PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &obj_index, dest );
  PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &align_row, dest );
  PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &align_col, dest );

  /* Add the contents to be added to the buffer */
  PLA_API_add_matrix_to_put_buffer( datatype, m, n, local_buf, lda, dest );

  /* Check if anything has come in for this node */
  PLA_API_Poll( );

  return PLA_SUCCESS;
}

/******************************************************************************/

static int put_sync_counter = 0;
static int request_sync_counter = 0;

int PLA_API_Poll ( )

/*--------------------------------------------------------------------------

Purpose :  Poll to see if messages have come in

----------------------------------------------------------------------------*/
{
  PLA_API_Poll_request( );
  PLA_API_Poll_put( );
}


int PLA_API_Poll_put ( )

/*--------------------------------------------------------------------------

Purpose :  Poll to see if messages have come in

----------------------------------------------------------------------------*/
{
  int 
    msgtype,
    m, n, lda,
    obj_index,
    align_row, align_col,
    local_ldim;

  PLA_Obj
    obj_dest = NULL;

  void
    *local_buf;

  MPI_Datatype
    datatype;

  if ( PLA_API_Accept_recv_put_buffer( ) ){
    /* New recv put_buffer has come in */
    while ( !PLA_API_End_of_recv_put_buffer ( ) ){
      /* Process next message */
      PLA_API_read_data_from_put_buffer( sizeof( int ), ( char * ) &msgtype );

      switch ( msgtype ){
      case PLA_API_TYPE_PUT:
	PLA_API_read_data_from_put_buffer( sizeof( int ), ( char * ) &m );
	PLA_API_read_data_from_put_buffer( sizeof( int ), ( char * ) &n );
	PLA_API_read_data_from_put_buffer( sizeof( int ), 
					   ( char * ) &obj_index );
	PLA_API_read_data_from_put_buffer( sizeof( int ), 
					   ( char * ) &align_row );
	PLA_API_read_data_from_put_buffer( sizeof( int ), 
					   ( char * ) &align_col );

	/* Take view into appropriate object */
	PLA_Obj_view( open_objects_list[ obj_index ],
		      m, n, align_row, align_col, &obj_dest );

	/* Add the contents of the message to the local contents
           of the view */
	PLA_API_read_contents_from_buffer_to_obj( m, n, obj_dest );

	break;
      case PLA_API_TYPE_SYNC:
	put_sync_counter ++;

	break;
      case PLA_API_TYPE_GET:
	/* Enter put destination information in the buffer */
	/*	PLA_API_read_data_from_put_buffer( sizeof( MPI_Datatype ), 
		( char *) &datatype ); */
	PLA_API_read_data_from_put_buffer( sizeof( int ), 
					    ( char * ) &obj_index );
	PLA_API_read_data_from_put_buffer( sizeof( int ), ( char * ) &m );
	PLA_API_read_data_from_put_buffer( sizeof( int ), ( char * ) &n );
	PLA_API_read_data_from_put_buffer( sizeof( void * ),
					    ( char * ) &local_buf );
	PLA_API_read_data_from_put_buffer( sizeof( int ), ( char * ) &lda );

	PLA_Obj_datatype( open_objects_list[ obj_index ], &datatype );

	/* Add the contents of the object to the buffer */
	PLA_API_add_matrix_from_recv_buffer_to_local( datatype,
               m, n, local_buf, lda );

	break;
      default:
	PLA_Abort( "API illegal type", __LINE__, __FILE__ );
      }
    }
  }

  PLA_Obj_free( &obj_dest );

  return PLA_SUCCESS;
}

/******************************************************************************/
int PLA_API_Poll_request ( )

/*--------------------------------------------------------------------------

Purpose :  Poll to see if messages have come in

----------------------------------------------------------------------------*/
{
  int 
    msgtype,
    m, n, lda,
    obj_index, source,
    align_row, align_col,
    typesize,
    local_ldim,
    temp;

  PLA_Obj
    obj_source = NULL;

  void
    *local_buf;

  MPI_Datatype 
    datatype;

  if ( PLA_API_Accept_recv_request_buffer( ) ){
    /* New recv put_buffer has come in */
    while ( !PLA_API_End_of_recv_request_buffer ( ) ){
      /* Process next message */
      PLA_API_read_data_from_request_buffer( sizeof( int ), 
					     ( char * ) &msgtype );
      switch ( msgtype ){
      case PLA_API_TYPE_PUT:
	PLA_Abort( "API illegal type", __LINE__, __FILE__ );

      case PLA_API_TYPE_SYNC:
	request_sync_counter ++;
	break;

      case PLA_API_TYPE_GET:
	PLA_API_read_data_from_request_buffer( sizeof( int ), ( char * ) &m );
	PLA_API_read_data_from_request_buffer( sizeof( int ), ( char * ) &n );
	PLA_API_read_data_from_request_buffer( sizeof( int ), 
					       ( char * ) &obj_index );
	PLA_API_read_data_from_request_buffer( sizeof( int ), 
					       ( char * ) &align_row );
	PLA_API_read_data_from_request_buffer( sizeof( int ), 
					       ( char * ) &align_col );
	PLA_API_read_data_from_request_buffer( sizeof( void * ),
						  ( char * ) &local_buf );
	PLA_API_read_data_from_request_buffer( sizeof( int ), 
					       ( char * )  &lda );
	PLA_API_read_data_from_request_buffer( sizeof( int ), 
					       ( char * )  &source );

	/* Take view into appropriate object */
	PLA_Obj_view( open_objects_list[ obj_index ],
		      m, n, align_row, align_col, &obj_source );

	/* Create room to send current message */
	PLA_Obj_datatype( obj_source, &datatype );
	MPI_Type_size( datatype, &typesize );
	PLA_API_put_buffer_make_room( 5 * sizeof( int ) + 
				      sizeof( void * ) +
				      m * n * typesize, source );
			   
	/* Enter put destination information in the buffer */
	temp = PLA_API_TYPE_GET;
	PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &temp,
				        source );
	/*	PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &datatype,
		source ); */
	PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &obj_index, 
				        source );
	PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &m, 
				        source );
	PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &n,
				        source );
	PLA_API_add_data_to_put_buffer( sizeof( void * ), 
				        ( char * ) &local_buf, source );
	PLA_API_add_data_to_put_buffer( sizeof( int ),  ( char * ) &lda,
				        source );

	/* Add the contents of the object to the buffer */
	PLA_API_add_contents_from_obj_to_put_buffer( obj_source, source );

	break;

      default:
	PLA_Abort( "API illegal type", __LINE__, __FILE__ );
      }
    }
  }

  PLA_Obj_free( &obj_source );

  return PLA_SUCCESS;
}

/******************************************************************************/

int pla2_api_sync_put( int repost )

/*--------------------------------------------------------------------------

Purpose : Synchronizes all nodes w.r.t. the given object and completes
all pending "API_axpy" operations on the object.

----------------------------------------------------------------------------*/
{
  int 
    dest, 
    nprocs,
    me, 
    ready,
    temp;

  MPI_Comm
    comm = MPI_COMM_NULL;

  MPI_Status status;

  PLA_Base_comm_1d( &comm );
  MPI_Comm_size ( comm, &nprocs );
  MPI_Comm_rank ( comm, &me );

  for ( dest = 0; dest < nprocs; dest++ ){
    if ( dest == me ) continue;
    /* Make sure there is room for synchronization signal to 
       be added to buffer */
    PLA_API_put_buffer_make_room( sizeof( int ), dest );

    /* Add synchronization signal to send buffer */
    temp = PLA_API_TYPE_SYNC;
    PLA_API_add_data_to_put_buffer( sizeof( int ), ( char * ) &temp, dest );

    /* Flush the asynchronous send buffer */
    PLA_API_put_buffer_flush( dest, repost );
  }

  /* Update synchronization counter and wait for it to reach zero */
  put_sync_counter -= nprocs-1;

  {
#if CHECK_FOR_STALL == 1
    double time, total_time;

    total_time = 0.0;
    time = MPI_Wtime();
#endif

    while ( put_sync_counter < 0 ) {
      PLA_API_Poll( );

#if CHECK_FOR_STALL == 1
      if ( ( MPI_Wtime() - time ) > STALL_INTERVAL ){
	total_time += MPI_Wtime() - time;
	printf("%d stalled in pla2_api_sync_put for %lf sec.\n", me, total_time );
	fflush( stdout );

	time = MPI_Wtime();
      }
#endif
    }
  }

  /* Just to make sure everything is in sync! */
#if CHECK_FOR_STALL == 1
  {
    int me;
    MPI_Comm_rank( comm, &me );

    printf("%d calling barrier in pla2_api_sync_put\n", me );
    fflush( stdout );
  }
#endif

  MPI_Barrier( comm );

#if CHECK_FOR_STALL == 1
  {
    int me;
    MPI_Comm_rank( comm, &me );

    printf("%d passed barrier in pla2_api_sync_put\n", me );
    fflush( stdout );
  }
#endif

  return PLA_SUCCESS;
}

/******************************************************************************/

int pla2_api_sync_request( int repost )

/*--------------------------------------------------------------------------

Purpose : Synchronizes the request buffers of 
    all nodes w.r.t. the given object and completes
    all pending "API_axpy" operations on the object.

----------------------------------------------------------------------------*/
{
  int 
    dest, 
    nprocs,
    me, 
    ready,
    temp;

  MPI_Comm
    comm = MPI_COMM_NULL;

  MPI_Status status;

  PLA_Base_comm_1d( &comm );
  MPI_Comm_size ( comm, &nprocs );
  MPI_Comm_rank ( comm, &me );

  for ( dest = 0; dest < nprocs; dest++ ){
    if ( dest == me ) continue;
    /* Make sure there is room for synchronization signal to 
       be added to buffer */
    PLA_API_request_buffer_make_room( sizeof( int ), dest );

    /* Add synchronization signal to send buffer */
    temp = PLA_API_TYPE_SYNC;
    PLA_API_add_data_to_request_buffer( sizeof( int ), ( char * ) &temp, 
				        dest );

    /* Flush the asynchronous send buffer */
    PLA_API_request_buffer_flush( dest, repost );
  }

  /* Update synchronization counter and wait for it to reach zero */
  request_sync_counter -= nprocs-1;

  {
#if CHECK_FOR_STALL == 1
    double time, total_time;

    total_time = 0.0;
    time = MPI_Wtime();
#endif
    while ( request_sync_counter < 0 ) {
      PLA_API_Poll( );

#if CHECK_FOR_STALL == 1
      if ( ( MPI_Wtime() - time ) > STALL_INTERVAL ){
	total_time += MPI_Wtime() - time;
	printf("%d stalled in pla2_api_sync_request for %lf sec.\n", 
	       me, total_time );
	fflush( stdout );
	time = MPI_Wtime();
      }
#endif
    }
  }

  /* Just to make sure everything is in sync! */
#if CHECK_FOR_STALL == 1
  {
    int me;
    MPI_Comm_rank( comm, &me );

    printf("%d calling barrier in pla2_api_sync_request", me );
    fflush( stdout );
  }
#endif

  MPI_Barrier( comm );

#if CHECK_FOR_STALL == 1
  {
    int me;
    MPI_Comm_rank( comm, &me );

    printf("%d passed barrier in pla2_api_sync_request\n", me );
    fflush( stdout );
  }
#endif

  return PLA_SUCCESS;
}


/******************************************************************************/

int PLA_API_put_buffer_flush( int dest, int repost )

/*--------------------------------------------------------------------------

Purpose : Flush the asynchronous put buffer to node dest

IN       dest      destination of buffer to be flushed

----------------------------------------------------------------------------*/
{
  int 
    ready;

  MPI_Comm
    comm = MPI_COMM_NULL;

  MPI_Status status;

  MPI_Request request;

  PLA_Base_comm_1d( &comm );

  if ( put_buffer_list[ dest ].cur_length != 0 ) {
    ready = FALSE;
    
    /* Wait until destination is ready for message.  
       Poll for incoming messages while doing this */
    while ( TRUE ){
#if CHECK_FOR_STALL == 1
      double time, total_time;
      int me;
      
      MPI_Comm_rank( comm, &me );

      total_time = 0.0;
      time = MPI_Wtime();
#endif

      MPI_Test( &put_buffer_list[ dest ].request, &ready, &status );
      if ( ready ) break;
      PLA_API_Poll_put( );

#if CHECK_FOR_STALL == 1
      if ( ( MPI_Wtime() - time ) > STALL_INTERVAL ){
	total_time += MPI_Wtime() - time;
	printf("%d stalled in PLA_API_put_buffer_flush for %lf sec.\n", 
	       me, total_time );
	fflush( stdout );

	time = MPI_Wtime();
      }
#endif
    }

    /* Create request for next message indicating the
       destination is ready to receive */
    if ( repost )
      MPI_Irecv( &dummy, 1, MPI_INT, dest, PUT_READY_TAG, comm,
                 &put_buffer_list[ dest ].request );

    /* Send the asynchronous send buffer to dest. */
#if CHECK_FOR_STALL == 1
    {
      int me;

      MPI_Comm_rank( comm, &me );
      printf("%d sending put buffer to %d\n", me, dest );
      fflush( stdout );
    }
#endif

    MPI_Isend( put_buffer_list[ dest ].buffer, 
              put_buffer_list[ dest ].cur_length,
	      MPI_CHAR, dest, PUT_TAG, comm, &request );

    while ( TRUE ){
      MPI_Test( &request, &ready, &status );
      if ( ready ) break;
      PLA_API_Poll_put( );
    }

#if CHECK_FOR_STALL == 1
    {
      int me;

      MPI_Comm_rank( comm, &me );
      printf("%d back from sending put buffer to %d\n", me, dest );
      fflush( stdout );
    }
#endif

    /* reset length of asynchronous send buffer */
    put_buffer_list[ dest ].cur_length = 0;
  }

  return PLA_SUCCESS;
}


int PLA_API_request_buffer_flush( int dest, int repost )

/*--------------------------------------------------------------------------

Purpose : Flush the asynchronous request buffer to node dest

IN       dest      destination of buffer to be flushed

----------------------------------------------------------------------------*/
{
  int 
    ready;

  MPI_Comm
    comm = MPI_COMM_NULL;

  MPI_Status 
    status;

  MPI_Request
    request;

  PLA_Base_comm_1d( &comm );

  if ( request_buffer_list[ dest ].cur_length != 0 ) {
    ready = FALSE;
    
    /* Wait until destination is ready for message.  
       Poll for incoming messages while doing this */
    while ( TRUE ){
#if CHECK_FOR_STALL == 1
      double time, total_time;
      int me;

      MPI_Comm_rank( comm, &me );

      total_time = 0.0;
      time = MPI_Wtime();
#endif

      MPI_Test( &request_buffer_list[ dest ].request, &ready, &status );
      if ( ready ) break;
      PLA_API_Poll( );

#if CHECK_FOR_STALL == 1
      if ( ( MPI_Wtime() - time ) > STALL_INTERVAL ){
	total_time += MPI_Wtime() - time;
	printf("%d stalled in pla2_api_sync_put for %lf sec.\n", 
	       me, total_time );
	fflush( stdout );
	time = MPI_Wtime();
      }
#endif
    }

    /* Create request for next message indicating the
       destination is ready to receive */
    if ( repost )
      MPI_Irecv( &dummy, 1, MPI_INT, dest, REQUEST_READY_TAG, comm,
                 &request_buffer_list[ dest ].request );

    /* Send the asynchronous send buffer to dest. */
#if CHECK_FOR_STALL == 1
    {
      int me;

      MPI_Comm_rank( comm, &me );
      printf("%d sending request buffer to %d\n", me, dest );
      fflush( stdout );
    }
#endif
    MPI_Isend( request_buffer_list[ dest ].buffer, 
              request_buffer_list[ dest ].cur_length,
	      MPI_CHAR, dest, REQUEST_TAG, comm, &request );

    while ( TRUE ){
      MPI_Test( &request, &ready, &status );
      if ( ready ) break;
      PLA_API_Poll( );
    }

#if CHECK_FOR_STALL == 1
    {
      int me;

      MPI_Comm_rank( comm, &me );
      printf("%d back from sending put buffer to %d\n", me, dest );
      fflush( stdout );
    }
#endif

    /* reset length of asynchronous send buffer */
    request_buffer_list[ dest ].cur_length = 0;
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_put_buffer_make_room( int size, int dest )

/*--------------------------------------------------------------------------

Purpose : Check if room for size more bites in asynchronous send
          buffer to node dest.  If not, flush.

IN       size      size to be added
IN       dest      destination of buffer

----------------------------------------------------------------------------*/
{
  size += SIZE_OF_UNIT;

  if ( size > MAX_PUT_BUFFER_SIZE )
    PLA_Abort( "size too large in PLA_API_put_buffer_make_room", __LINE__, __FILE__ );

  if ( put_buffer_list[ dest ].cur_length + size > MAX_PUT_BUFFER_SIZE ) {
    PLA_API_put_buffer_flush( dest, TRUE );
    put_buffer_list[ dest ].cur_length = 0;
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_request_buffer_make_room( int size, int dest )

/*--------------------------------------------------------------------------

Purpose : Check if room for size more bites in asynchronous send
          buffer to node dest.  If not, flush.

IN       size      size to be added
IN       dest      destination of buffer

----------------------------------------------------------------------------*/
{
  size += SIZE_OF_UNIT;

  if ( size > MAX_REQUEST_BUFFER_SIZE )
    PLA_Abort( "size too large in PLA_API_request_buffer_make_room", __LINE__, __FILE__) ;

  if ( request_buffer_list[ dest ].cur_length + size > 
                MAX_REQUEST_BUFFER_SIZE ) {
    PLA_API_request_buffer_flush( dest, TRUE );
    request_buffer_list[ dest ].cur_length = 0;
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_add_data_to_put_buffer( int size, char *data, int dest )

/*--------------------------------------------------------------------------

Purpose : add integer to end of asynchronous send buffer.

IN       data      data to be added
IN       dest      destination of buffer

----------------------------------------------------------------------------*/
{
  int 
    i;
  char *tempp;

  tempp = ( char * ) put_buffer_list[ dest ].buffer + 
                     put_buffer_list[ dest ].cur_length;

  for ( i=0; i<size; i++ )
    *tempp++ = *data++;

  put_buffer_list[ dest ].cur_length += size;

/*
  for ( i=0; i<size; i++ ){
    * ( ( char * ) ( put_buffer_list[ dest ].buffer 
      + put_buffer_list[ dest ].cur_length ) ) = *data++;
    put_buffer_list[ dest ].cur_length++;
  }
*/
  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_add_data_to_request_buffer( int size, char *data, int dest )

/*--------------------------------------------------------------------------

Purpose : add integer to end of asynchronous send buffer.

IN       data      data to be added
IN       dest      destination of buffer

----------------------------------------------------------------------------*/
{
  int 
    i;
  char *tempp;

  tempp = ( char * ) request_buffer_list[ dest ].buffer + 
                     request_buffer_list[ dest ].cur_length;

  for ( i=0; i<size; i++ )
    *tempp++ = *data++;

  request_buffer_list[ dest ].cur_length += size;

/*  for ( i=0; i<size; i++ ){
    * ( ( char * ) ( request_buffer_list[ dest ].buffer 
      + request_buffer_list[ dest ].cur_length ) ) = *data++;
    request_buffer_list[ dest ].cur_length++;
  } */

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_add_matrix_to_put_buffer( 
     MPI_Datatype datatype, int m, int n, void *local_buf, int lda, int dest )

/*--------------------------------------------------------------------------

Purpose : add matrix to end of asynchronous send buffer.

----------------------------------------------------------------------------*/
{
  char 
    *tempp_from, *tempp_to;
  int 
    i, j,
    typesize, cur_length;

  /* align */
  cur_length = put_buffer_list[ dest ].cur_length;
  cur_length = ( ( cur_length % SIZE_OF_UNIT ) == 0 ? 
		 cur_length : ( cur_length / SIZE_OF_UNIT + 1 ) * SIZE_OF_UNIT );
  put_buffer_list[ dest ].cur_length = cur_length;

  tempp_to = put_buffer_list[ dest ].buffer + 
             put_buffer_list[ dest ].cur_length;

  MPI_Type_size( datatype, &typesize );

  put_buffer_list[ dest ].cur_length += typesize * m * n;

  if ( put_buffer_list[ dest ].cur_length > MAX_PUT_BUFFER_SIZE )
    PLA_Abort( "sending too much", __LINE__, __FILE__ );

  lda = lda * typesize;
  m = m * typesize;
  for ( j=0; j<n; j++ ){
    tempp_from = ( char * ) local_buf + j * lda;
    for ( i=0; i<m; i++ )
      *tempp_to++ = *tempp_from++;
  }

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_read_data_from_put_buffer( int size, char *data )

/*--------------------------------------------------------------------------

Purpose : read data from location in recv buffer.

IN       size      length of data in characters
IN       data      address of data

----------------------------------------------------------------------------*/
{
  int 
    i;
  char 
    *tempp;

  tempp = ( char * ) recv_put_buffer + recv_put_buffer_done;

  for ( i=0; i<size; i++ )
    *data++ = *tempp++;

  recv_put_buffer_done += size;

  return PLA_SUCCESS;
}

int PLA_API_read_data_from_request_buffer( int size, char *data )

/*--------------------------------------------------------------------------

Purpose : read integer from location in recv buffer.

IN       size      size of data
IN       data      address of data

----------------------------------------------------------------------------*/
{

  int 
    i;
  char
    *tempp;

  tempp = ( char * ) recv_request_buffer + recv_request_buffer_done;

  for ( i=0; i<size; i++ )
    *data++ = *tempp++;

  recv_request_buffer_done += size;

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_read_contents_from_buffer_to_obj( int m, int n, PLA_Obj obj )

/*--------------------------------------------------------------------------

Purpose : read matrix from current location in recv buffer.

----------------------------------------------------------------------------*/
{
  void 
    *local_buf;

  int 
    local_m, local_n, local_ldim, typesize, cur_length;

  MPI_Datatype
    datatype;

  char 
    *tempp;


  PLA_Obj_local_buffer( obj, ( void ** ) &local_buf );
  PLA_Obj_local_ldim( obj, &local_ldim );
  PLA_Obj_datatype( obj, &datatype );
  MPI_Type_size( datatype, &typesize );

  /* align */
  cur_length = recv_put_buffer_done;
  cur_length = ( ( cur_length % SIZE_OF_UNIT ) == 0 ? 
		 cur_length : ( cur_length / SIZE_OF_UNIT + 1 ) * SIZE_OF_UNIT );
  recv_put_buffer_done = cur_length;

  tempp = ( char * ) recv_put_buffer + recv_put_buffer_done;

  PLA_Obj_add_to_local_contents( PLA_NO_TRANS, m, n, 
              ( void * ) tempp,
              m, 1, obj );

  recv_put_buffer_done += m * n * typesize;

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_add_contents_from_obj_to_put_buffer( PLA_Obj obj, int dest )

/*--------------------------------------------------------------------------

Purpose : read matrix from object into put buffer

----------------------------------------------------------------------------*/
{
  char
    *bufp;

  int 
    local_m, local_n, dummy, typesize, cur_length; 

  MPI_Datatype
    datatype;

  cur_length = put_buffer_list[ dest ].cur_length;
  cur_length = ( ( cur_length % SIZE_OF_UNIT ) == 0 ? 
		 cur_length : ( cur_length / SIZE_OF_UNIT + 1 ) * SIZE_OF_UNIT );
  put_buffer_list[ dest ].cur_length = cur_length;

  bufp = ( char *) put_buffer_list[ dest ].buffer + 
                   put_buffer_list[ dest ].cur_length;

  PLA_Obj_local_length( obj, &local_m );
  PLA_Obj_local_width( obj, &local_n );

  
  PLA_Obj_get_local_contents( obj, PLA_NO_TRANSPOSE, &dummy, &dummy,
			      bufp, local_m, 1 );

  PLA_Obj_datatype( obj, &datatype );
  MPI_Type_size( datatype, &typesize );
	
  put_buffer_list[ dest ].cur_length += local_m * local_n * typesize;

  return PLA_SUCCESS;
}



int PLA_API_read_matrix_from_buffer( 
     MPI_Datatype datatype, int m, int n, void *local_buf, int lda )

/*--------------------------------------------------------------------------

Purpose : read matrix from current location in recv buffer.

----------------------------------------------------------------------------*/
{
  char 
    *tempp_from, *tempp_to;
  int 
    i, j,
    typesize, cur_length;

  /* align */
  cur_length = recv_put_buffer_done;
  cur_length = ( ( cur_length % SIZE_OF_UNIT ) == 0 ? 
		 cur_length : ( cur_length / SIZE_OF_UNIT + 1 ) * SIZE_OF_UNIT );
  recv_put_buffer_done = cur_length;

  tempp_from = ( char * ) recv_put_buffer + recv_put_buffer_done;

  MPI_Type_size( datatype, &typesize );

  lda = lda * typesize;
  m = m * typesize;
  for ( j=0; j<n; j++ ){
    tempp_to = ( char * ) local_buf + j * lda;
    for ( i=0; i<m; i++ )
      *tempp_to++ = *tempp_from++;
  }

  recv_put_buffer_done += m * n;

  return PLA_SUCCESS;
}

/******************************************************************************/

int PLA_API_Accept_recv_put_buffer( )

/*--------------------------------------------------------------------------

Purpose : Check if asynchronous put buffer has arrived, and accept it
          if it has.

----------------------------------------------------------------------------*/
{
  int 
    source,
    value;

  MPI_Comm
    comm = MPI_COMM_NULL;

  MPI_Status
    status;

  PLA_Base_comm_1d( &comm );

  MPI_Iprobe( MPI_ANY_SOURCE, PUT_TAG, comm, &value, &status );

  source = status.MPI_SOURCE;
    
  if ( value ){
#if CHECK_FOR_STALL == 1
    int me;

    MPI_Comm_rank( comm, &me );
    printf("%d received put buffer from %d\n", me, source );
    fflush( stdout );
#endif

    MPI_Recv( recv_put_buffer, MAX_PUT_BUFFER_SIZE, MPI_CHAR,
	      source, PUT_TAG, comm, &status );

    MPI_Get_count( &status, MPI_CHAR, &recv_put_buffer_length );

    recv_put_buffer_done = 0;

    MPI_Send( &dummy, 1, MPI_INT, 
	      source, PUT_READY_TAG, comm );
  }

  return value;
}


int PLA_API_Accept_recv_request_buffer( )

/*--------------------------------------------------------------------------

Purpose : Check if asynchronous request buffer has arrived, and accept it
          if it has.

----------------------------------------------------------------------------*/
{
  int 
    source,
    value;

  MPI_Comm
    comm = MPI_COMM_NULL;

  MPI_Status
    status;

  PLA_Base_comm_1d( &comm );

  MPI_Iprobe( MPI_ANY_SOURCE, REQUEST_TAG, comm, &value, &status );

  source = status.MPI_SOURCE;
    
  if ( value ){
#if CHECK_FOR_STALL == 1
    int me;

    MPI_Comm_rank( comm, &me );
    printf("%d received request buffer from %d\n", me, source );
    fflush( stdout );
#endif

    MPI_Recv( recv_request_buffer, MAX_REQUEST_BUFFER_SIZE, MPI_CHAR,
	      source, REQUEST_TAG, comm, &status );

    MPI_Get_count( &status, MPI_CHAR, &recv_request_buffer_length );

    recv_request_buffer_done = 0;

    MPI_Send( &dummy, 1, MPI_INT, 
	      source, REQUEST_READY_TAG, comm );
  }

  return value;
}

/******************************************************************************/

int PLA_API_End_of_recv_put_buffer ( ) 

/*--------------------------------------------------------------------------

Purpose : Check if end of receive buffer.

----------------------------------------------------------------------------*/
{
  return ( recv_put_buffer_done >= recv_put_buffer_length );
}

/******************************************************************************/

int PLA_API_End_of_recv_request_buffer ( ) 

/*--------------------------------------------------------------------------

Purpose : Check if end of receive buffer.

----------------------------------------------------------------------------*/
{
  return ( recv_request_buffer_done >= recv_request_buffer_length );
}



int PLA_Obj_index_in_open_objects( PLA_Obj obj, int *index )
{
  int 
    i, 
    count;

  count = 0;
  *index = -1;
  for ( i=0; i<MAX_OPEN_OBJS; i++ ){
    if ( obj == open_objects_list[ i ] ){
      *index = i;
      break;
    }
  }

  if ( *index == -1 )
    PLA_Abort( "error determining index of open object", __LINE__, __FILE__ );

  return PLA_SUCCESS;
}


/******************************************************************************/

int PLA_Obj_add_to_local_contents   (
              int      trans,                 int      rows_in_buf,
              int      cols_in_buf,           void     *buf,
              int      leading_dim_buf,       int      stride_buf,
              PLA_Obj obj )

/*----------------------------------------------------------------------------

Purpose : Add to the local data from the given object.


IN        trans             indicates whether to transpose data
IN        rows_in_buf       row dimension of data buffer
IN        cols_in_buf       column dimension of data buffer
IN        buf               address of data buffer
IN        leading_dim       leading dimension of data buffer
IN        stride_buf        stride of data buffer where data is put
IN/OUT    obj               global object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    m, n, i, j, ldim, stride;
  char 
    *buf_local, *buf_obj, *tempp_local, *tempp_obj;
  MPI_Datatype 
    datatype;
  
  PLA_Obj_local_length( obj, &m );
  PLA_Obj_local_width ( obj, &n );
  PLA_Obj_local_ldim  ( obj, &ldim );    
  PLA_Obj_datatype    ( obj, &datatype );

  if ( datatype == MPI_DOUBLE ){
    double *buf_local, *buf_obj, *tempp_local, *tempp_obj;
    
    buf_local = (double *) buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

    if ( trans == PLA_NO_TRANS ){
      for ( j=0; j<n; j++ ){
	tempp_local = buf_local + j*leading_dim_buf;
	tempp_obj   = buf_obj   + j*ldim;
	for ( i=0; i<m; i++ )
	  *tempp_obj++ += *tempp_local++;
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = buf_local + j*leading_dim_buf;
	tempp_obj   = buf_obj + j;
	for ( i=0; i<n; i++ ) {
	  *tempp_obj += *tempp_local;
	  tempp_obj += ldim;
	}
      }
    }
  }
  else if ( datatype == MPI_FLOAT) {
    float *buf_local, *buf_obj, *tempp_local, *tempp_obj;
    
    buf_local = (float *) buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

    if ( trans == PLA_NO_TRANS ){
      for ( j=0; j<n; j++ ){
	tempp_local = buf_local + j*leading_dim_buf;
	tempp_obj   = buf_obj   + j*ldim;
	for ( i=0; i<m; i++ )
	  *tempp_obj++ += *tempp_local++;
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = buf_local + j*leading_dim_buf;
	tempp_obj   = buf_obj + j;
	for ( i=0; i<n; i++ ) {
	  *tempp_obj += *tempp_local;
	  tempp_obj += ldim;
	}
      }
    }
  }
  else if ( datatype == MPI_COMPLEX) {
    float *buf_local, *buf_obj, *tempp_local, *tempp_obj;
    
    buf_local = (float *) buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

    if ( trans == PLA_NO_TRANS ){
      for ( j=0; j<n; j++ ){
	tempp_local = buf_local + j*leading_dim_buf*2;
	tempp_obj   = buf_obj   + j*ldim*2;
	for ( i=0; i<m*2; i++ )
	  *tempp_obj++ += *tempp_local++;
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = buf_local + j*leading_dim_buf*2;
	tempp_obj   = buf_obj + j*2;
	for ( i=0; i<n; i++ ) {
	  *tempp_obj++ += *tempp_local++;
	  *tempp_obj++ += *tempp_local++;
	  tempp_obj += (ldim-1)*2;
	}
      }
    }
  }
  else if ( datatype == MPI_DOUBLE_COMPLEX) {
    double *buf_local, *buf_obj, *tempp_local, *tempp_obj;
    
    buf_local = (double *) buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

    if ( trans == PLA_NO_TRANS ){
      for ( j=0; j<n; j++ ){
	tempp_local = buf_local + j*leading_dim_buf*2;
	tempp_obj   = buf_obj   + j*ldim*2;
	for ( i=0; i<m*2; i++ )
	  *tempp_obj++ += *tempp_local++;
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = buf_local + j*leading_dim_buf*2;
	tempp_obj   = buf_obj + j*2;
	for ( i=0; i<n; i++ ) {
	  *tempp_obj++ += *tempp_local++;
	  *tempp_obj++ += *tempp_local++;
	  tempp_obj += (ldim-1)*2;
	}
      }
    }
  }
  else {
    PLA_Abort( "PLA_Obj_add_to_local_contents: datatype not implemented", __LINE__, __FILE__ );
  }

  return value;
}

/******************************************************************************/

int PLA_Obj_add_from_local_contents   (
              int      trans,                 int      rows_in_buf,
              int      cols_in_buf,           void     *buf,
              int      leading_dim_buf,       int      stride_buf,
              PLA_Obj obj )

/*----------------------------------------------------------------------------

Purpose : Add the local object contents to the local buffer


IN        trans             indicates whether to transpose data
IN        rows_in_buf       row dimension of data buffer
IN        cols_in_buf       column dimension of data buffer
IN        buf               address of data buffer
IN        leading_dim       leading dimension of data buffer
IN        stride_buf        stride of data buffer where data is put
IN/OUT    obj               global object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    m, n, i, j, ldim, stride;
  char 
    *buf_local, *buf_obj, *tempp_local, *tempp_obj;
  MPI_Datatype 
    datatype;
  
  PLA_Obj_local_length( obj, &m );
  PLA_Obj_local_width ( obj, &n );
  PLA_Obj_local_ldim  ( obj, &ldim );    
  PLA_Obj_datatype    ( obj, &datatype );

  if ( datatype == MPI_DOUBLE) {
    double *buf_local, *buf_obj, *tempp_local, *tempp_obj;
    
    buf_local = (double *) buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

    if ( trans == PLA_NO_TRANS ){
      for ( j=0; j<n; j++ ){
	tempp_local = buf_local + j*leading_dim_buf;
	tempp_obj   = buf_obj   + j*ldim;
	for ( i=0; i<m; i++ )
	  *tempp_local++ += *tempp_obj++;
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = buf_local + j*leading_dim_buf;
	tempp_obj   = buf_obj + j;
	for ( i=0; i<n; i++ ) {
	  *tempp_local++ += *tempp_obj;
	  tempp_obj += ldim;
	}
      }
    }
  }
  else if ( datatype == MPI_FLOAT) {
    float *buf_local, *buf_obj, *tempp_local, *tempp_obj;
    
    buf_local = (float *) buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

    if ( trans == PLA_NO_TRANS ){
      for ( j=0; j<n; j++ ){
	tempp_local = buf_local + j*leading_dim_buf;
	tempp_obj   = buf_obj   + j*ldim;
	for ( i=0; i<m; i++ )
	  *tempp_local++ += *tempp_obj++;
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = buf_local + j*leading_dim_buf;
	tempp_obj   = buf_obj + j;
	for ( i=0; i<n; i++ ) {
	  *tempp_local++ += *tempp_obj;
	  tempp_obj += ldim;
	}
      }
    }
  }
  else if ( datatype == MPI_DOUBLE_COMPLEX) {
    double *buf_local, *buf_obj, *tempp_local, *tempp_obj;
    
    buf_local = ( double *) buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

    if ( trans == PLA_NO_TRANS ){
      for ( j=0; j<n; j++ ){
	tempp_local = buf_local + j*leading_dim_buf*2;
	tempp_obj   = buf_obj   + j*ldim*2;
	for ( i=0; i<2*m; i++ )
	  *tempp_local++ += *tempp_obj++;
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = buf_local + j*leading_dim_buf*2;
	tempp_obj   = buf_obj + j*2;
	for ( i=0; i<n; i++ ) {
	  *tempp_local++ += *tempp_obj++;
	  *tempp_local++ += *tempp_obj++;
	  tempp_obj += (ldim-1)*2;
	}
      }
    }
  }
  else if ( datatype == MPI_COMPLEX) {
    float *buf_local, *buf_obj, *tempp_local, *tempp_obj;
    
    buf_local = ( float *) buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

    if ( trans == PLA_NO_TRANS ){
      for ( j=0; j<n; j++ ){
	tempp_local = buf_local + j*leading_dim_buf*2;
	tempp_obj   = buf_obj   + j*ldim*2;
	for ( i=0; i<2*m; i++ )
	  *tempp_local++ += *tempp_obj++;
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = buf_local + j*leading_dim_buf*2;
	tempp_obj   = buf_obj + j*2;
	for ( i=0; i<n; i++ ) {
	  *tempp_local++ += *tempp_obj++;
	  *tempp_local++ += *tempp_obj++;
	  tempp_obj += (ldim-1)*2;
	}
      }
    }
  }
  else {
    PLA_Abort( "PLA_Obj_add_from_local_contents: datatype not implemented", __LINE__, __FILE__ );
  }

  return value;
}

int PLA_API_add_matrix_from_recv_buffer_to_local( 
	       MPI_Datatype datatype, int m, int n, 
               void *local_buf, int lda )
{
  int i, j, typesize, cur_length;

  /* align */
  cur_length = recv_put_buffer_done;
  cur_length = ( ( cur_length % SIZE_OF_UNIT ) == 0 ? 
		 cur_length : ( cur_length / SIZE_OF_UNIT + 1 ) * SIZE_OF_UNIT );
  recv_put_buffer_done = cur_length;

  MPI_Type_size( datatype, &typesize );

  if ( datatype == MPI_DOUBLE ) {
    double *buf_local, *buf_recv, *tempp_local, *tempp_recv;
    
    buf_local = (double *) local_buf;
    buf_recv = ( double * ) ( ( char * ) recv_put_buffer + recv_put_buffer_done );

    for ( j=0; j<n; j++ ){
      tempp_local = buf_local + j*lda;
      tempp_recv  = buf_recv   + j*m;
      for ( i=0; i<m; i++ )
	*tempp_local++ += *tempp_recv++;
    }

    recv_put_buffer_done += m * n * typesize;
  }
  else if ( datatype == MPI_FLOAT ) {
    float *buf_local, *buf_recv, *tempp_local, *tempp_recv;
    
    buf_local = (float *) local_buf;
    buf_recv = ( float * ) ( ( char * ) recv_put_buffer + recv_put_buffer_done );

    for ( j=0; j<n; j++ ){
      tempp_local = buf_local + j*lda;
      tempp_recv  = buf_recv   + j*m;
      for ( i=0; i<m; i++ )
	*tempp_local++ += *tempp_recv++;
    }

    recv_put_buffer_done += m * n * typesize;
  }
  else if ( datatype == MPI_DOUBLE_COMPLEX ) {
    double *buf_local, *buf_recv, *tempp_local, *tempp_recv;
    
    buf_local = (double *) local_buf;
    buf_recv = ( double * ) ( ( char * ) recv_put_buffer + recv_put_buffer_done );

    for ( j=0; j<n; j++ ){
      tempp_local = buf_local + j*lda*2;
      tempp_recv  = buf_recv   + j*m*2;
      for ( i=0; i<2*m; i++ )
	*tempp_local++ += *tempp_recv++;
    }

    recv_put_buffer_done += m * n * typesize;
  }
  else if ( datatype == MPI_COMPLEX ) {
    float *buf_local, *buf_recv, *tempp_local, *tempp_recv;
    
    buf_local = (float *) local_buf;
    buf_recv = ( float * ) ( ( char * ) recv_put_buffer + recv_put_buffer_done );

    for ( j=0; j<n; j++ ){
      tempp_local = buf_local + j*lda*2;
      tempp_recv  = buf_recv   + j*m*2;
      for ( i=0; i<2*m; i++ )
	*tempp_local++ += *tempp_recv++;
    }

    recv_put_buffer_done += m * n * typesize;
  }
  else {
    PLA_Abort( "datatype not yet implemented", __LINE__, __FILE__ );
  }
}
    

int PLA_API_add_obj_to_open_objects_list( PLA_Obj obj )
{
  int i;

  if ( number_open_objects >= MAX_OPEN_OBJS )
    PLA_Abort( "maximum number of open objects exceeded", __LINE__, __FILE__ );

  /* Find an unused entry in the open_objects_list */
  for ( i=0; i<MAX_OPEN_OBJS; i++ )
    if ( open_objects_list[ i ] == NULL ) break;

  if ( i == MAX_OPEN_OBJS )
    PLA_Abort( "maximum number of open objects exceeded", __LINE__, __FILE__ );

  open_objects_list[ i ] = obj;

  number_open_objects++;
  
  return PLA_SUCCESS;
}  


int PLA_API_delete_obj_from_open_objects_list( PLA_Obj obj )
{
  int i;

  /* Free the entry in the open objects list */
  PLA_Obj_index_in_open_objects( obj, &i );
  open_objects_list[ i ] = NULL;

  number_open_objects--;
  return PLA_SUCCESS;
}  


int PLA_API_number_open_objects( int *num_open )
{
  *num_open = number_open_objects;
  
  return PLA_SUCCESS;
}

int PLA_API_flush_sync_messages( )
{
  int
    dummy, i, nprocs;

  MPI_Comm
    comm = MPI_COMM_NULL;

  MPI_Status
    status;

  PLA_Base_comm_1d( &comm );

  MPI_Comm_size( comm, &nprocs );

  for ( i=0; i<nprocs-1; i++ )
    MPI_Recv( &dummy, 1, MPI_INT,
	      MPI_ANY_SOURCE, 999, comm, &status );
  for ( i=0; i<nprocs-1; i++ )
    MPI_Recv( &dummy, 1, MPI_INT,
	      MPI_ANY_SOURCE, 997, comm, &status );
}
