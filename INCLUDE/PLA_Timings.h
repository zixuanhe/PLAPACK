/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/*---------------------------------------------------
 Structures and defines for timing data
---------------------------------------------------*/

/* Local BLAS1 Operations */
#define PLA_LOCAL_ASUM_TIMING  10
#define PLA_LOCAL_AXPY_TIMING  11
#define PLA_LOCAL_COPY_TIMING  12
#define PLA_LOCAL_DOT_TIMING   13
#define PLA_LOCAL_IAMAX_TIMING 14
#define PLA_LOCAL_NRM2_TIMING  15
#define PLA_LOCAL_SCAL_TIMING  16
#define PLA_LOCAL_SWAP_TIMING  17

/* Local BLAS2 Operations */
#define PLA_LOCAL_GEMV_TIMING  20
#define PLA_LOCAL_GER_TIMING   21
#define PLA_LOCAL_HEMV_TIMING  22
#define PLA_LOCAL_SYMV_TIMING  23
#define PLA_LOCAL_SYR_TIMING   24
#define PLA_LOCAL_TRMV_TIMING  25
#define PLA_LOCAL_TRSV_TIMING  26

/* Local BLAS3 Operations */
#define PLA_LOCAL_GEMM_TIMING  30
#define PLA_LOCAL_HERK_TIMING  31
#define PLA_LOCAL_SYMM_TIMING  32
#define PLA_LOCAL_SYRK_TIMING  33
#define PLA_LOCAL_TRMM_TIMING  34
#define PLA_LOCAL_TRSM_TIMING  35
#define PLA_LOCAL_GEMM_TIMING_IMBALANCE  36

/* BLAS1 Operations */
#define PLA_ASUM_TIMING        40
#define PLA_AXPY_TIMING        41

#define PLA_DOT_TIMING         43
#define PLA_IAMAX_TIMING       44
#define PLA_NRM2_TIMING        45
#define PLA_SCAL_TIMING        46
#define PLA_SWAP_TIMING        47

/* Communications Operations */
#define PLA_COPY_TIMING  50

/* MPI Operations */
#define PLA_MPI_SEND_TIMING    60
#define PLA_MPI_RECV_TIMING    61
#define PLA_MPI_IRECV_TIMING   62
#define PLA_MPI_WAIT_TIMING    63
#define PLA_MPI_BCAST_TIMING   64


#define PLA_LOCAL_LU_TIMING 90

#define PLA_TIMINGS_SIZE  100

extern double * PLA_TIMINGS;

/* Timing Routines */

int PLA_Timings_initialize( );

int PLA_Timings_print( );

int PLA_Timings_average( );

