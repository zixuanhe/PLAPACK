/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/* 
   This implementation of the Application Interface keeps a list
   of open objects, the open_objects_list.
*/

/******************************************************************************/

int PLA_API_Init_open_objects_list( );

/*--------------------------------------------------------------------------

Purpose : Initialize open objects list

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_Init_put_buffers( );

/*--------------------------------------------------------------------------

Purpose : Initialize put buffers

----------------------------------------------------------------------------*/

/************************************************************************/

int PLA_API_Init_request_buffers( );

/*--------------------------------------------------------------------------

Purpose : Initialize request buffers

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_Free_open_objects_list( );

/*--------------------------------------------------------------------------

Purpose : Free open objects list

----------------------------------------------------------------------------*/

/******************************************************************************/

int   PLA_API_Free_put_buffers( );

/*--------------------------------------------------------------------------

Purpose : Free put buffers

----------------------------------------------------------------------------*/

/******************************************************************************/

int   PLA_API_Free_request_buffers( );

/*--------------------------------------------------------------------------

Purpose : Free request buffers

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_request( MPI_Datatype, int, int, void *, int, int, int, int, int );

/*--------------------------------------------------------------------------

Purpose : enter a request for data on the appropriate request buffer 
          to be sent to processor "dest"

----------------------------------------------------------------------------*/

/******************************************************************************/
int PLA_API_put  ( MPI_Datatype, int, int, void *, int, int, int, int, int );

/*--------------------------------------------------------------------------

Purpose : enter data in the appropriate put buffer to be sent to 
          processor "dest"

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_Poll ( );

/*--------------------------------------------------------------------------

Purpose :  Poll to see if messages have come in

----------------------------------------------------------------------------*/

int PLA_API_Poll_put ( );

/*--------------------------------------------------------------------------

Purpose :  Poll to see if messages have come in

----------------------------------------------------------------------------*/

/******************************************************************************/
int PLA_API_Poll_request ( );

/*--------------------------------------------------------------------------

Purpose :  Poll to see if messages have come in

----------------------------------------------------------------------------*/

/******************************************************************************/

int pla2_api_sync_put( int );

/*--------------------------------------------------------------------------

Purpose : Synchronizes all nodes w.r.t. the given object and completes
all pending "API_axpy" operations on the object.

----------------------------------------------------------------------------*/

/******************************************************************************/

int pla2_api_sync_request( int );

/*--------------------------------------------------------------------------

Purpose : Synchronizes the request buffers of 
    all nodes w.r.t. the given object and completes
    all pending "API_axpy" operations on the object.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_put_buffer_flush( int, int );

/*--------------------------------------------------------------------------

Purpose : Flush the asynchronous put buffer to node dest

IN       dest      destination of buffer to be flushed

----------------------------------------------------------------------------*/

int PLA_API_request_buffer_flush( int, int );

/*--------------------------------------------------------------------------

Purpose : Flush the asynchronous request buffer to node dest

IN       dest      destination of buffer to be flushed

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_put_buffer_make_room( int, int );

/*--------------------------------------------------------------------------

Purpose : Check if room for size more bites in asynchronous send
          buffer to node dest.  If not, flush.

IN       size      size to be added
IN       dest      destination of buffer

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_request_buffer_make_room( int, int );

/*--------------------------------------------------------------------------

Purpose : Check if room for size more bites in asynchronous send
          buffer to node dest.  If not, flush.

IN       size      size to be added
IN       dest      destination of buffer

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_add_data_to_put_buffer( int, char *, int );

/*--------------------------------------------------------------------------

Purpose : add integer to end of asynchronous send buffer.

IN       data      data to be added
IN       dest      destination of buffer

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_add_data_to_request_buffer( int, char *, int );

/*--------------------------------------------------------------------------

Purpose : add integer to end of asynchronous send buffer.

IN       data      data to be added
IN       dest      destination of buffer

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_add_matrix_to_put_buffer( 
                 MPI_Datatype, int, int, void *, int, int );

/*--------------------------------------------------------------------------

Purpose : add matrix to end of asynchronous send buffer.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_read_data_from_put_buffer( int, char * );

/*--------------------------------------------------------------------------

Purpose : read data from location in recv buffer.

IN       size      length of data in characters
IN       data      address of data

----------------------------------------------------------------------------*/

int PLA_API_read_data_from_request_buffer( int, char * );

/*--------------------------------------------------------------------------

Purpose : read integer from location in recv buffer.

IN       size      size of data
IN       data      address of data

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_read_pntr_from_request_buffer( void ** );

/*--------------------------------------------------------------------------

Purpose : read integer from location in recv buffer.

IN       data      address of pointer

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_read_pntr_from_put_buffer( void ** );

/*--------------------------------------------------------------------------

Purpose : read integer from location in recv buffer.

IN       data      address of pointer

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_read_contents_from_buffer_to_obj( int, int, PLA_Obj );

/*--------------------------------------------------------------------------

Purpose : read matrix from current location in recv buffer.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_add_contents_from_obj_to_put_buffer( PLA_Obj, int );

/*--------------------------------------------------------------------------

Purpose : read matrix from object into put buffer

----------------------------------------------------------------------------*/


int PLA_API_read_matrix_from_buffer( 
     MPI_Datatype, int, int, void *, int );

/*--------------------------------------------------------------------------

Purpose : read matrix from current location in recv buffer.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_Accept_recv_put_buffer( );

/*--------------------------------------------------------------------------

Purpose : Check if asynchronous put buffer has arrived, and accept it
          if it has.

----------------------------------------------------------------------------*/


int PLA_API_Accept_recv_request_buffer( );

/*--------------------------------------------------------------------------

Purpose : Check if asynchronous request buffer has arrived, and accept it
          if it has.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_End_of_recv_put_buffer ( );

/*--------------------------------------------------------------------------

Purpose : Check if end of receive buffer.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_End_of_recv_request_buffer ( );

/*--------------------------------------------------------------------------

Purpose : Check if end of receive buffer.

----------------------------------------------------------------------------*/

int PLA_Obj_index_in_open_objects( PLA_Obj obj, int *index );


/******************************************************************************/

int PLA_Obj_add_to_local_contents ( int, int, int, void *, int, int, PLA_Obj );

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

/******************************************************************************/

int PLA_Obj_add_from_local_contents(int, int, int, void *, int, int, PLA_Obj );

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

int PLA_API_add_matrix_from_recv_buffer_to_local( 
                 MPI_Datatype, int, int, void *, int );
    

int PLA_API_add_obj_to_open_objects_list( PLA_Obj );

int PLA_API_delete_obj_from_open_objects_list( PLA_Obj );

int PLA_API_number_open_objects( int * );

