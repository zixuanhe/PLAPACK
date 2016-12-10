/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

static MPI_Comm
  base_comm    = MPI_COMM_NULL,
  base_comm_1d = MPI_COMM_NULL;

static int    
  PLA_initialized = FALSE;

int  PLA_Comm_1D_to_2D (MPI_Comm  comm1d,
                        int       rows,
                        int       cols,
                        MPI_Comm  *comm2d)
{
int
   dims[2],
   periods[2]={TRUE, TRUE},
   reorder=FALSE,
   error=PLA_SUCCESS;

   dims[0] = rows; dims[1] = cols;          
   error+=MPI_Cart_create(comm1d, 2, dims, periods, reorder, comm2d);

   return (error);
}


int  PLA_Comm_1D_to_2D_ratio (MPI_Comm  comm1d,
                              double    ratio,
                              MPI_Comm  *comm2d)
/* 
   Want to find rows and colums such that 
              rows/columns = ratio
   or at least a close approxximation.  Noting that 
   columns=num_proc/rows and rearranging
              rows*rows  = ratio*num_proc
   So find rows such that 
         min |rows*rows - ratio*num_proc|
   and (num_proc%rows)=0
*/
{
int
   error = 0,
   num_proc,
   target_rows,
   minimum,
   temp,
   factor,
   rows;

   error+= MPI_Comm_size  (comm1d, &num_proc);

   rows    = num_proc;
   minimum = num_proc*num_proc;
   target_rows = ratio*(float)num_proc + .5; 
   for ( factor=1; factor<=num_proc ; factor++ ) {
      if ( num_proc%factor == 0) { /* this is an integer factor */
         temp = factor*factor - target_rows;
         if (temp<0) temp = -temp;
         if (temp<minimum) {
	     minimum = temp;
             rows    = factor;
         }
      }
   }
   error += PLA_Comm_1D_to_2D (comm1d,
                               rows,
                               num_proc/rows,
                               comm2d);

    return(error);
}


int  PLA_Init ( MPI_Comm  comm )
{
  if ( PLA_initialized ) 
    PLA_Abort( "PLA_Init called more than once", __LINE__, __FILE__ );

  PLA_initialized = TRUE;

  if ( MPI_COMM_NULL == comm ) 
    PLA_Abort( "Null communicator", __LINE__, __FILE__);

  MPI_Comm_dup ( comm, &base_comm );
  
  pla_comm_create_1d_from_2d( base_comm, &base_comm_1d );

#ifdef PLA_MPI_COMPLEX
    MPI_Type_contiguous (2, MPI_FLOAT, &MPI_COMPLEX);
    MPI_Type_commit (&MPI_COMPLEX);
#endif

#ifdef PLA_MPI_DOUBLE_COMPLEX
    MPI_Type_contiguous (2, MPI_DOUBLE, &MPI_DOUBLE_COMPLEX);
    MPI_Type_commit (&MPI_DOUBLE_COMPLEX);
#endif

  
  return PLA_SUCCESS;
}

int PLA_Initialized ( int *initialized ) 
{
  return PLA_initialized;
}

int PLA_Base_comm ( MPI_Comm *comm )
{
   *comm = base_comm;

   return PLA_SUCCESS;
}

int PLA_Base_comm_1d ( MPI_Comm *comm )
{
   *comm = base_comm_1d;

   return PLA_SUCCESS;
}

int   PLA_Finalize ()
{
  if ( !PLA_initialized )
    PLA_Abort( "PLAPACK not yet initialized", __LINE__, __FILE__ );

  PLA_initialized = FALSE;

  MPI_Comm_free ( &base_comm );
  MPI_Comm_free ( &base_comm_1d );

  base_comm    = MPI_COMM_NULL;
  base_comm_1d = MPI_COMM_NULL;

#ifdef PLA_MPI_COMPLEX
    MPI_Type_free (&MPI_COMPLEX);
#endif

#ifdef PLA_MPI_COMPLEX
    MPI_Type_free (&MPI_DOUBLE_COMPLEX);
#endif

  return PLA_SUCCESS;
}

int pla_comm_create_1d_from_2d( MPI_Comm  comm2d,
				MPI_Comm  *comm1d )
{
int   
  error=0,
  my_cvect_index,
  dims[2], 
  periods[2], 
  coords[2];
  /*
     Create a 1D communicator.  This must be in column major order and
     created from the cartisian grid passed in the calling sequence.
     They will be ordered based on the number my_cvect_index.  So we
     just have to make sure it is increasing in column major order.
  */

  error = MPI_Cart_get(comm2d, 2, dims, periods, coords);
  if (error) 
    PLA_Abort( "Error returned from MPI_Cart_get", __LINE__, __FILE__ );

  my_cvect_index = coords[1] * dims[0] + coords[0];
  error = MPI_Comm_split( comm2d, 0, my_cvect_index, comm1d);
  if (error) 
    PLA_Abort( "Error returned from MPI_Comm_split", __LINE__, __FILE__ );

  return PLA_SUCCESS;
}
