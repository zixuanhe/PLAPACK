/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"


#ifndef min
#define min(x,y) (x <= y ? x : y)
#endif


int PLA_Temp_create (int nb, int zero_or_one, PLA_Template *templ_ret)

{
  int      
    dims[2],
    remain_dims[2],
    error;

  MPI_Comm     
    comm;

  PLA_Template
    templ;

  dims[0]=dims[1]=0;

  /* 
     Allocate Memory for the platemplate.  
  */

  templ = ( PLA_Template ) 
    PLA_malloc( (size_t) sizeof( struct PLA_Template_struct ) );

  if ( templ == NULL)
    PLA_Abort( "malloc failed", __LINE__, __FILE__ );

  /* 
     Get the global communicator.
  */

  PLA_Base_comm( &comm );
  if ( MPI_COMM_NULL==comm ) 
    PLA_Abort( "PLA_Init not called before PLA_Temp_create", __LINE__, __FILE__ );
  
  /* 
     Store the parameters passed.
  */

  templ->nb       = nb;
  templ->comm_all = comm;
  if ( zero_or_one != 0 && zero_or_one != 1 )
    PLA_Abort( "Error in PLA_temp_create: zero_or_one", __LINE__, __FILE__ );

  templ->zero_or_one = zero_or_one;

  /* 
     Set the other 2 communicators, row, column, and 1D.
     remain_dims tells which dimensions to keep and which to throw
     away.  For instance remain_dims[0]=1, remain_dims[1]=0 keeps the 
     rows and throws away the columns creating a col communicator.
  */
  
  remain_dims[0]=1;
  remain_dims[1]=0;
  templ->comm_col = MPI_COMM_NULL;
  error = MPI_Cart_sub ( comm, remain_dims, &(templ->comm_col)); 
  if ( error )
    PLA_Abort( "Error creating column communicator", __LINE__, __FILE__ );

  remain_dims[0]=0;
  remain_dims[1]=1;
  templ->comm_row = MPI_COMM_NULL;
  error = MPI_Cart_sub( comm, remain_dims, &( templ->comm_row));
  if ( error )
    PLA_Abort( "Error creating row communicator", __LINE__, __FILE__ );

  /* 
     To compute the row dimension, look at the size of the column 
     communicator.  To compute the size of the column dimentin
     look at the size of the row communicator.
  */
  error = MPI_Comm_size ( templ->comm_col, &(templ->nprows) );
  if ( error )
    PLA_Abort( "Error computing nprows", __LINE__, __FILE__ );

  error = MPI_Comm_size ( templ->comm_row, &(templ->npcols) );
  if ( error )
    PLA_Abort( "Error computing npcols", __LINE__, __FILE__ );

  /* 
     Store the position of this node within the processing grid.
     The postion in the row communicator is the column index
     and vice versa.
  */
  error = MPI_Comm_rank ( templ->comm_col, &(templ->myrow) );
  if ( error )
    PLA_Abort( "Error computing myrow", __LINE__, __FILE__ );

  error = MPI_Comm_rank ( templ->comm_row, &(templ->mycol) );
  if ( error )
    PLA_Abort( "Error computing mycol", __LINE__, __FILE__ );

  error = PLA_Base_comm_1d( &(templ->comm_all) );
  if ( error )
    PLA_Abort( "Error returned from PLA_Base_comm_1d", __LINE__, __FILE__ );

  error = MPI_Comm_rank( templ->comm_all, &(templ->me) );
  if ( error )
    PLA_Abort( "Error computing me", __LINE__, __FILE__ );
  
  /* 
     Now set the number of times this platemplate has been used.  
     This number will be incremented when a matrix is allocated
     using this platemplate.  It will be decremented when that
     matrix is released.  The platemplate will not be released 
     until this number has been decremented to zero.

  */
  templ->times_used = 0;

  *templ_ret = templ;

  return PLA_SUCCESS;
}



int  PLA_Temp_free( PLA_Template *templ )

{
  /* 
     I guess check for NULL just to be safe.
  */
  if ( NULL != *templ ){
    if ( (*templ)->times_used == 0) {
      MPI_Comm_free ( &( (*templ)->comm_col ) );
      MPI_Comm_free ( &( (*templ)->comm_row ) );
      
      PLA_free ( *templ );
    }
    *templ = NULL;
  }

  return PLA_SUCCESS;
}

int PLA_Temp_comm_col_size( PLA_Template templ, int *nprows )
{
  *nprows = ( templ->nprows );

  return PLA_SUCCESS;
}

int PLA_Temp_comm_row_size( PLA_Template templ, int *npcols )
{
  *npcols = ( templ->npcols );

  return PLA_SUCCESS;
}

int PLA_Temp_comm_all_size( PLA_Template templ, int *nprocs )
{
  *nprocs = ( templ->nprows ) * ( templ->npcols );

  return PLA_SUCCESS;
}

int PLA_Temp_nb( PLA_Template templ, int *nb )
{
  *nb = ( templ->nb );

  return PLA_SUCCESS;
}

int PLA_Temp_comm_all_info( PLA_Template templ,
			     MPI_Comm    *comm,
			     int         *me,
			     int         *nprocs )
{
   int return_value = 0;

   return_value += PLA_Temp_comm_all( templ, comm);
   return_value += PLA_Temp_comm_all_rank( templ, me);
   return_value += PLA_Temp_comm_all_size( templ, nprocs );

   return return_value;
}

int PLA_Temp_comm_all( PLA_Template templ,
                       MPI_Comm     *comm )
{
   *comm = templ->comm_all;
   
   return PLA_SUCCESS;
}

int PLA_Temp_comm_all_rank( PLA_Template  templ,
			     int           *me )
{
  *me = templ->me;

  return PLA_SUCCESS;
}

int PLA_Temp_comm_row_info( PLA_Template   templ,
			     MPI_Comm       *comm_row,
			     int            *mycol,
			     int            *npcols )
{
  *comm_row = templ->comm_row;
  *mycol    = templ->mycol;
  *npcols   = templ->npcols;

  return PLA_SUCCESS;
}

int PLA_Temp_comm_col( PLA_Template  templ,
			MPI_Comm      *comm_col)
{
  *comm_col = templ->comm_col;

  return PLA_SUCCESS;
}

int PLA_Temp_comm_col_info( PLA_Template  templ,
			     MPI_Comm      *comm_col,
			     int           *myrow,
			     int           *nprows )
{
  *comm_col = templ->comm_col;
  *myrow    = templ->myrow;
  *nprows   = templ->nprows;

  return PLA_SUCCESS;
}

int PLA_Temp_comm_col_rank( PLA_Template templ,
			     int          *myrow )
{
  *myrow = templ->myrow;

  return PLA_SUCCESS;
}

int PLA_Temp_comm_row( PLA_Template  templ,
			MPI_Comm      *comm_row)
{
  *comm_row = templ->comm_row;

  return PLA_SUCCESS;
}

int PLA_Temp_comm_row_rank( PLA_Template templ,
			     int          *mycol )
{
  *mycol = templ->mycol;

  return PLA_SUCCESS;
}

int PLA_Temp_zero_or_one( PLA_Template templ,
			   int *zero_or_one)
{
  *zero_or_one = ( templ->zero_or_one );

  return PLA_SUCCESS;
}

