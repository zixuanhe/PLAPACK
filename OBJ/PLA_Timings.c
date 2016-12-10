/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

double * PLA_TIMINGS = NULL;

/***************************************************************************/

int PLA_Timings_initialize( )
     
/*----------------------------------------------------------------------------

Purpose : Initialize structures for detailed timing

----------------------------------------------------------------------------*/
{
  int i;

  if (PLA_TIMINGS == NULL) {
    PLA_TIMINGS = (double *) malloc( sizeof(double) * PLA_TIMINGS_SIZE );
  }

  for (i = 0; i<PLA_TIMINGS_SIZE; i++)
    PLA_TIMINGS[i] = 0;

  PLA_ERROR_CHECKING = TRUE;

  return PLA_SUCCESS;
}

/***************************************************************************/

int PLA_Timings_print()
     
/*----------------------------------------------------------------------------

Purpose : Print timing data

----------------------------------------------------------------------------*/
{
  int i;

  double total_time = 0.0;

  if (PLA_TIMINGS == NULL)
    return PLA_SUCCESS;

  /* Local BLAS1 Timings */
  if ( PLA_TIMINGS[PLA_LOCAL_ASUM_TIMING] != 0 )
    printf("Local_asum  time = %f\n", PLA_TIMINGS[PLA_LOCAL_ASUM_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_AXPY_TIMING] != 0 )
    printf("Local_axpy  time = %f\n", PLA_TIMINGS[PLA_LOCAL_AXPY_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_COPY_TIMING] != 0 )
    printf("Local_copy  time = %f\n", PLA_TIMINGS[PLA_LOCAL_COPY_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_DOT_TIMING] != 0 )
    printf("Local_dot   time = %f\n", PLA_TIMINGS[PLA_LOCAL_DOT_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_IAMAX_TIMING] != 0 )
    printf("Local_iamax time = %f\n", PLA_TIMINGS[PLA_LOCAL_IAMAX_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_SCAL_TIMING] != 0 )
    printf("Local_scal  time = %f\n", PLA_TIMINGS[PLA_LOCAL_SCAL_TIMING]);

  /* Local BLAS2 Timings */
  if ( PLA_TIMINGS[PLA_LOCAL_GEMV_TIMING] != 0 )
    printf("Local_Gemv  time = %f\n", PLA_TIMINGS[PLA_LOCAL_GEMV_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_GER_TIMING] != 0 )
    printf("Local_Ger   time = %f\n", PLA_TIMINGS[PLA_LOCAL_GER_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_HEMV_TIMING] != 0 )
    printf("Local_Hemv  time = %f\n", PLA_TIMINGS[PLA_LOCAL_HEMV_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_SYMV_TIMING] != 0 )
    printf("Local_Symv  time = %f\n", PLA_TIMINGS[PLA_LOCAL_SYMV_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_SYR_TIMING] != 0 )
    printf("Local_Syr   time = %f\n", PLA_TIMINGS[PLA_LOCAL_SYR_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_TRMV_TIMING] != 0 )
    printf("Local_Trmv  time = %f\n", PLA_TIMINGS[PLA_LOCAL_TRMV_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_TRSV_TIMING] != 0 )
    printf("Local_Trsv  time = %f\n", PLA_TIMINGS[PLA_LOCAL_TRSV_TIMING]);

  /* Local BLAS3 Timings */
  if ( PLA_TIMINGS[PLA_LOCAL_GEMM_TIMING] != 0 )
    printf("Local_Gemm  time = %f\n", PLA_TIMINGS[PLA_LOCAL_GEMM_TIMING]);
  if ( PLA_TIMINGS[PLA_LOCAL_GEMM_TIMING_IMBALANCE] != 0 )
    printf("Local_Gemm_imb  time = %f\n", PLA_TIMINGS[PLA_LOCAL_GEMM_TIMING_IMBALANCE]);
  if ( PLA_TIMINGS[PLA_LOCAL_TRSM_TIMING] != 0 )
    printf("Local_Trsm  time = %f\n", PLA_TIMINGS[PLA_LOCAL_TRSM_TIMING]);

  /* MPI Timings */
  if ( PLA_TIMINGS[PLA_MPI_SEND_TIMING] != 0 )
    printf("MPI_Send    time = %f\n", PLA_TIMINGS[PLA_MPI_SEND_TIMING]);
  if ( PLA_TIMINGS[PLA_MPI_RECV_TIMING] != 0 )
    printf("MPI_Recv    time = %f\n", PLA_TIMINGS[PLA_MPI_RECV_TIMING]);
  if ( PLA_TIMINGS[PLA_MPI_IRECV_TIMING] != 0 )
    printf("MPI_Irecv   time = %f\n", PLA_TIMINGS[PLA_MPI_IRECV_TIMING]);
  if ( PLA_TIMINGS[PLA_MPI_WAIT_TIMING] != 0 )
    printf("MPI_Wait    time = %f\n", PLA_TIMINGS[PLA_MPI_WAIT_TIMING]);
  if ( PLA_TIMINGS[PLA_MPI_BCAST_TIMING] != 0 )
    printf("MPI_Bcast   time = %f\n", PLA_TIMINGS[PLA_MPI_BCAST_TIMING]);

  if ( PLA_TIMINGS[PLA_IAMAX_TIMING] != 0 )
    printf("Iamax time       = %f\n", PLA_TIMINGS[PLA_IAMAX_TIMING]);

  if ( PLA_TIMINGS[PLA_COPY_TIMING] != 0 )
    printf("Copy time        = %f\n", PLA_TIMINGS[PLA_COPY_TIMING]);

  if ( PLA_TIMINGS[PLA_LOCAL_LU_TIMING] != 0 )
    printf("Local_lu    time = %f\n", PLA_TIMINGS[PLA_LOCAL_LU_TIMING]);

  /* Print total time captured */

  for (i = 0; i<PLA_TIMINGS_SIZE; i++)
    total_time += PLA_TIMINGS[i];

  printf("**TOTAL**   time = %f\n", total_time);

  return PLA_SUCCESS;
}

/***************************************************************************/

int PLA_Timings_average()
     
/*----------------------------------------------------------------------------

Purpose : Compute averages across all nodes

----------------------------------------------------------------------------*/
{
  int nprocs, i;

  if (PLA_TIMINGS == NULL)
    return PLA_SUCCESS;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  MPI_Reduce( PLA_TIMINGS, PLA_TIMINGS, PLA_TIMINGS_SIZE, MPI_DOUBLE, MPI_SUM, 
                  /* root = */ 0, MPI_COMM_WORLD);

  for (i = 0; i<PLA_TIMINGS_SIZE; i++)
    PLA_TIMINGS[i] = PLA_TIMINGS[i] / nprocs;
   
  return PLA_SUCCESS;
}

