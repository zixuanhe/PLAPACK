/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#if MANUFACTURE == CRAY
#include <string.h>
#include <fortran.h>
#endif

/*
   FORTRAN strings are handled specially on CRAY systems.  All
   other systems should use the first definition.
 */

#if MANUFACTURE == CRAY
  #define PLA_CP2FCD( name ) _cptofcd ( name, strlen ( name ) )
#else  
  #define PLA_CP2FCD( name ) name
#endif

/*
   The calling convention from C to Fortran differs depending on platform.
   The only two options we have encountered are appending an underscore
   and doing nothing.
 */
#if     MACHINE_TYPE == PARAGON || MANUFACTURE == SUN || MANUFACTURE == PC || MANUFACTURE == SGI
  #define PLA_FOR2C( name ) name ## _
#elif MANUFACTURE == CRAY || MACHINE_TYPE == SP2 || MANUFACTURE == IBM || MANUFACTURE == HP
  #define PLA_FOR2C( name ) name
#endif


/* 
   The following definitions stem from the fact that on CRAY platforms 
   floats are 64 bits.  So even when one is using MPI_DOUBLE in a C
   program or MPI_DOUBLE_COMPLEX in a fortran program, one makes calls
   to the single precision BLAS, e.g. SGEMM ( ... ).  

   The only two categories here are MANUFACTURE (sic) == CRAY and != CRAY.
 */


/*
 * Level 1 BLAS
 */


#if MANUFACTURE != CRAY

  /*
   * Single precision
   */

  #define PLA_srotg( a, b, c, s )                         PLA_FOR2C ( srotg  ) ( a, b, c, s )
  #define PLA_srot(   n,        x, incx, y, incy, c, s )  PLA_FOR2C ( srot   ) ( n,        x, incx, y, incy, c, s )
  #define PLA_sswap(  n,        x, incx, y, incy )        PLA_FOR2C ( sswap  ) ( n,        x, incx, y, incy )
  #define PLA_scopy(  n,        x, incx, y, incy )        PLA_FOR2C ( scopy  ) ( n,        x, incx, y, incy )
  #define PLA_saxpy(  n, alpha, x, incx, y, incy )        PLA_FOR2C ( saxpy  ) ( n, alpha, x, incx, y, incy )
  #define PLA_sscal(  n, alpha, x, incx)                  PLA_FOR2C ( sscal  ) ( n, alpha, x, incx)
  #define PLA_sdot(   n,        x, incx, y, incy )        PLA_FOR2C ( sdot   ) ( n,        x, incx, y, incy )
  #define PLA_snrm2(  n,        x, incx )                 PLA_FOR2C ( snrm2  ) ( n,        x, incx )
  #define PLA_sasum(  n,        x, incx )                 PLA_FOR2C ( sasum  ) ( n,        x, incx )
  #define PLA_isamax( n,        x, incx )                 PLA_FOR2C ( isamax ) ( n,        x, incx )

  #define PLA_slacpy( u, m, n, a, lda, b, ldb )           PLA_FOR2C ( slacpy ) ( u, m, n, a, lda, b, ldb)

  /* prototypes */

#if     MACHINE_TYPE == PARAGON || MANUFACTURE == SUN || MANUFACTURE == PC || MANUFACTURE == SGI
  void  srotg_  ( float *, float *, float *, float * );
  void  srot_   ( int *, float *, int *, float *, int *, float *, float * );
  void  sswap_  ( int *, float *, int *, float *, int * );
  void  scopy_  ( int *, float *, int *, float *, int * );
  void  saxpy_  ( int *, float *, float *, int *, float *, int * );
  void  sscal_  ( int *, float *, float *, int * );
  float sdot_   ( int *, float *, int *, float *, int * );
  float snrm2_  ( int *, float *, int * );
  float sasum_  ( int *, float *, int * );
  int   isamax_ ( int *, float *, int * );
#elif MANUFACTURE == CRAY || MACHINE_TYPE == SP2 || MANUFACTURE == IBM || MANUFACTURE == HP
  void  srotg ( float *, float *, float *, float * );
  void  srot  ( int *, float *, int *, float *, int *, float *, float * );
  void  sswap ( int *, float *, int *, float *, int * );
  void  scopy ( int *, float *, int *, float *, int * );
  void  saxpy ( int *, float *, float *, int *, float *, int * );
  void  sscal ( int *, float *, float *, int * );
  float sdot  ( int *, float *, int *, float *, int * );
  float snrm2 ( int *, float *, int * );
  float sasum ( int *, float *, int * );
  int   isamax( int *, float *, int * );
#endif


  /*
   * Double precision
   */

  #define PLA_drotg( a, b, c, s )                         PLA_FOR2C ( drotg  ) ( a, b, c, s )
  #define PLA_drot(   n,        x, incx, y, incy, c, s )  PLA_FOR2C ( drot   ) ( n,        x, incx, y, incy, c, s )
  #define PLA_dswap(  n,        x, incx, y, incy )        PLA_FOR2C ( dswap  ) ( n,        x, incx, y, incy )
  #define PLA_dcopy(  n,        x, incx, y, incy )        PLA_FOR2C ( dcopy  ) ( n,        x, incx, y, incy )
  #define PLA_daxpy(  n, alpha, x, incx, y, incy )        PLA_FOR2C ( daxpy  ) ( n, alpha, x, incx, y, incy )
  #define PLA_dscal(  n, alpha, x, incx)                  PLA_FOR2C ( dscal  ) ( n, alpha, x, incx)
  #define PLA_ddot(   n,        x, incx, y, incy )        PLA_FOR2C ( ddot   ) ( n,        x, incx, y, incy )
  #define PLA_dnrm2(  n,        x, incx )                 PLA_FOR2C ( dnrm2  ) ( n,        x, incx )
  #define PLA_dasum(  n,        x, incx )                 PLA_FOR2C ( dasum  ) ( n,        x, incx )
  #define PLA_idamax( n,        x, incx )                 PLA_FOR2C ( idamax ) ( n,        x, incx )
 
  #define PLA_dlacpy( u, m, n, a, lda, b, ldb )           PLA_FOR2C ( dlacpy ) ( u, m, n, a, lda, b, ldb)

  /* prototype */ 

#if     MACHINE_TYPE == PARAGON || MANUFACTURE == SUN || MANUFACTURE == PC || MANUFACTURE == SGI
  void  drotg_  ( double *, double *, double *, double * );
  void  drot_   ( int *, double *, int *, double *, int *, double *, double * );
  void  dswap_  ( int *, double *, int *, double *, int * );
  void  dcopy_  ( int *, double *, int *, double *, int * );
  void  daxpy_  ( int *, double *, double *, int *, double *, int * );
  void  dscal_  ( int *, double *, double *, int * );
  double ddot_   ( int *, double *, int *, double *, int * );
  double dnrm2_  ( int *, double *, int * );
  double dasum_  ( int *, double *, int * );
  int   idamax_ ( int *, double *, int * );
#elif MANUFACTURE == CRAY || MACHINE_TYPE == SP2 || MANUFACTURE == IBM || MANUFACTURE == HP
  void  drotg ( double *, double *, double *, double * );
  void  drot  ( int *, double *, int *, double *, int *, double *, double * );
  void  dswap ( int *, double *, int *, double *, int * );
  void  dcopy ( int *, double *, int *, double *, int * );
  void  daxpy ( int *, double *, double *, int *, double *, int * );
  void  dscal ( int *, double *, double *, int * );
  double ddot  ( int *, double *, int *, double *, int * );
  double dnrm2 ( int *, double *, int * );
  double dasum ( int *, double *, int * );
  int   idamax( int *, double *, int * );
#endif

  /*
   * Single precision complex
   */
  #define PLA_cswap(  n,        x, incx, y, incy )        PLA_FOR2C ( cswap  ) ( n,        x, incx, y, incy )
  #define PLA_ccopy(  n,        x, incx, y, incy )        PLA_FOR2C ( ccopy  ) ( n,        x, incx, y, incy )
  #define PLA_caxpy(  n, alpha, x, incx, y, incy )        PLA_FOR2C ( caxpy  ) ( n, alpha, x, incx, y, incy )
  #define PLA_cscal(  n, alpha, x, incx)                  PLA_FOR2C ( cscal  ) ( n, alpha, x, incx)
  #define PLA_cdot(   n,        x, incx, y, incy )        PLA_FOR2C ( cdot   ) ( n,        x, incx, y, incy )
  #define PLA_cdotu(  n,        x, incx, y, incy )        PLA_FOR2C ( cdotu  ) ( n,        x, incx, y, incy )
  #define PLA_cdotc(  n,        x, incx, y, incy )        PLA_FOR2C ( cdotc  ) ( n,        x, incx, y, incy )
  #define PLA_scnrm2(  n,        x, incx )                 PLA_FOR2C ( scnrm2  ) ( n,        x, incx )
  #define PLA_scasum(  n,        x, incx )                 PLA_FOR2C ( scasum  ) ( n,        x, incx )
  #define PLA_icamax( n,        x, incx )                 PLA_FOR2C ( icamax ) ( n,        x, incx )

  #define PLA_clacpy( u, m, n, a, lda, b, ldb )           PLA_FOR2C ( clacpy ) ( u, m, n, a, lda, b, ldb)

  /* prototypes */

#if     MACHINE_TYPE == PARAGON || MANUFACTURE == SUN || MANUFACTURE == PC || MANUFACTURE == SGI
  void  crotg_ ( PLA_COMPLEX *, PLA_COMPLEX *, PLA_COMPLEX *, PLA_COMPLEX * );
  void  crot_ ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, PLA_COMPLEX * );
  void  cswap_( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
  void  ccopy_( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
  void  caxpy_( int *, PLA_COMPLEX *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
  void  cscal_( int *, PLA_COMPLEX *, PLA_COMPLEX *, int * );
  PLA_COMPLEX cdot_ ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
  float scnrm2_( int *, PLA_COMPLEX *, int * );
  float scasum_( int *, PLA_COMPLEX *, int * );
  int   icamax_ ( int *, PLA_COMPLEX *, int * );
#elif MANUFACTURE == CRAY || MACHINE_TYPE == SP2 || MANUFACTURE == IBM || MANUFACTURE == HP
  void  crotg ( PLA_COMPLEX *, PLA_COMPLEX *, PLA_COMPLEX *, PLA_COMPLEX * );
  void  crot  ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, PLA_COMPLEX * );
  void  cswap ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
  void  ccopy ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
  void  caxpy ( int *, PLA_COMPLEX *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
  void  cscal ( int *, PLA_COMPLEX *, PLA_COMPLEX *, int * );
  PLA_COMPLEX cdot_  ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
  float scnrm2 ( int *, PLA_COMPLEX *, int * );
  float scasum ( int *, PLA_COMPLEX *, int * );
  int   icamax( int *, PLA_COMPLEX *, int * );
#endif

  /*
   * Double precision complex
   */
  #define PLA_zswap(  n,        x, incx, y, incy )        PLA_FOR2C ( zswap  ) ( n,        x, incx, y, incy )
  #define PLA_zcopy(  n,        x, incx, y, incy )        PLA_FOR2C ( zcopy  ) ( n,        x, incx, y, incy )
  #define PLA_zaxpy(  n, alpha, x, incx, y, incy )        PLA_FOR2C ( zaxpy  ) ( n, alpha, x, incx, y, incy )
  #define PLA_zscal(  n, alpha, x, incx)                  PLA_FOR2C ( zscal  ) ( n, alpha, x, incx)
  #define PLA_zdot(   n,        x, incx, y, incy )        PLA_FOR2C ( zdot   ) ( n,        x, incx, y, incy )
  #define PLA_zdotu(  n,        x, incx, y, incy )        PLA_FOR2C ( zdotu  ) ( n,        x, incx, y, incy )
  #define PLA_zdotc(  n,        x, incx, y, incy )        PLA_FOR2C ( zdotc  ) ( n,        x, incx, y, incy )
  #define PLA_dznrm2(  n,        x, incx )                 PLA_FOR2C ( dznrm2  ) ( n,        x, incx )
  #define PLA_dzasum(  n,        x, incx )                 PLA_FOR2C ( dzasum  ) ( n,        x, incx )
  #define PLA_izamax( n,        x, incx )                 PLA_FOR2C ( izamax ) ( n,        x, incx )

  #define PLA_zlacpy( u, m, n, a, lda, b, ldb )           PLA_FOR2C ( zlacpy ) ( u, m, n, a, lda, b, ldb)

  /* prototypes */

#if     MACHINE_TYPE == PARAGON || MANUFACTURE == SUN || MANUFACTURE == PC || MANUFACTURE == SGI
  void  zrotg_ ( PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX * );
  void  zrot_  ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX * );
  void  zswap_ ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  void  zcopy_ ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  void  zaxpy_ ( int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  void  zscal_ ( int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  PLA_DOUBLE_COMPLEX zdot_  ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  double dznrm2_ ( int *, PLA_DOUBLE_COMPLEX *, int * );
  double dzasum_ ( int *, PLA_DOUBLE_COMPLEX *, int * );
  int   izamax_( int *, PLA_DOUBLE_COMPLEX *, int * );
#elif MANUFACTURE == CRAY || MACHINE_TYPE == SP2 || MANUFACTURE == IBM || MANUFACTURE == HP
  void  zrotg ( PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX * );
  void  zrot  ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX * );
  void  zswap ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  void  zcopy ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  void  zaxpy ( int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  void  zscal ( int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  PLA_DOUBLE_COMPLEX zdot_  ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  double dznrm2 ( int *, PLA_DOUBLE_COMPLEX *, int * );
  double dzasum ( int *, PLA_DOUBLE_COMPLEX *, int * );
  int   izamax ( int *, PLA_DOUBLE_COMPLEX *, int * );
#endif

#else  /* MANUFACTURE == CRAY */

  /*
   * Single precision
   */

  #define PLA_srotg( a, b, c, s )                         PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_srot(   n,        x, incx, y, incy, c, s )  PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_sswap(  n,        x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_scopy(  n,        x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_saxpy(  n, alpha, x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_sscal(  n, alpha, x, incx)                  PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_sdot(   n,        x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_snrm2(  n,        x, incx )                 PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_sasum(  n,        x, incx )                 PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_isamax( n,        x, incx )                 PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_slacpy( u, m, n, a, lda, b, ldb )           PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  /*
   * Double precision
   */
  #define PLA_drotg( a, b, c, s )                         SROTG  ( a, b, c, s )
  #define PLA_drot(   n,        x, incx, y, incy, c, s )  SROT   ( n,        x, incx, y, incy, c, s )
  #define PLA_dswap(  n,        x, incx, y, incy )        SSWAP  ( n,        x, incx, y, incy )
  #define PLA_dcopy(  n,        x, incx, y, incy )        SCOPY  ( n,        x, incx, y, incy )
  #define PLA_daxpy(  n, alpha, x, incx, y, incy )        SAXPY  ( n, alpha, x, incx, y, incy )
  #define PLA_dscal(  n, alpha, x, incx)                  SSCAL  ( n, alpha, x, incx)
  #define PLA_ddot(   n,        x, incx, y, incy )        SDOT   ( n,        x, incx, y, incy )
  #define PLA_dnrm2(  n,        x, incx )                 SNRM2  ( n,        x, incx )
  #define PLA_dasum(  n,        x, incx )                 SASUM  ( n,        x, incx )
  #define PLA_idamax( n,        x, incx )                 ISAMAX ( n,        x, incx )

  #define PLA_dlacpy( u, m, n, a, lda, b, ldb )           SLACPY ( PLA_CP2FCD ( u ), m, n, a, lda, b, ldb)
   
  /* prototypes */

  void  SROTG  ( double *, double *, double *, double * );
  void  SROT   ( int *, double *, int *, double *, int *, double *, double * );
  void  SSWAP  ( int *, double *, int *, double *, int * );
  void  SCOPY  ( int *, double *, int *, double *, int * );
  void  SAXPY  ( int *, double *, double *, int *, double *, int * );
  void  SSCAL  ( int *, double *, double *, int * );
  double SDOT   ( int *, double *, int *, double *, int * );
  double SNRM2  ( int *, double *, int * );
  double SASUM  ( int *, double *, int * );
  int   ISAMAX ( int *, double *, int * );

  /*
   * Single precision complex
   */
  #define PLA_cswap(  n,        x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_ccopy(  n,        x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_caxpy(  n, alpha, x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_cscal(  n, alpha, x, incx)                  PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_cdot(   n,        x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_cdotu(  n,        x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_cdotc(  n,        x, incx, y, incy )        PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_scnrm2(  n,        x, incx )                 PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_scasum(  n,        x, incx )                 PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_icamax( n,        x, incx )                 PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_clacpy( u, m, n, a, lda, b, ldb )           PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  /*
   * Double precision complex
   */
  #define PLA_zswap(  n,        x, incx, y, incy )        CSWAP  ( n,        x, incx, y, incy )
  #define PLA_zcopy(  n,        x, incx, y, incy )        CCOPY  ( n,        x, incx, y, incy )
  #define PLA_zaxpy(  n, alpha, x, incx, y, incy )        CAXPY  ( n, alpha, x, incx, y, incy )
  #define PLA_zscal(  n, alpha, x, incx)                  CSCAL  ( n, alpha, x, incx)
  #define PLA_zdot(   n,        x, incx, y, incy )        CDOT   ( n,        x, incx, y, incy )
  #define PLA_zdotu(  n,        x, incx, y, incy )        CDOTU  ( n,        x, incx, y, incy )
  #define PLA_zdotc(  n,        x, incx, y, incy )        CDOTC  ( n,        x, incx, y, incy )
  #define PLA_dznrm2(  n,        x, incx )                SCNRM2  ( n,        x, incx )
  #define PLA_dzasum(  n,        x, incx )                 SCASUM  ( n,        x, incx )
  #define PLA_izamax( n,        x, incx )                 ICAMAX ( n,        x, incx )

  #define PLA_zlacpy( u, m, n, a, lda, b, ldb )           CLACPY ( PLA_CP2FCD ( u ), m, n, a, lda, b, ldb)

  /* prototypes */

  void  CROTG  ( PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, 
		 PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX * );
  void  CROT   ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, 
		 PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX * );
  void  CSWAP  ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  void  CCOPY  ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  void  CAXPY  ( int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, 
                                                PLA_DOUBLE_COMPLEX *, int * );
  void  CSCAL  ( int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  PLA_DOUBLE_COMPLEX CDOT   ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  double SCNRM2  ( int *, PLA_DOUBLE_COMPLEX *, int * );
  double SCASUM  ( int *, PLA_DOUBLE_COMPLEX *, int * );
  int    ICAMAX ( int *, PLA_DOUBLE_COMPLEX *, int * );

#endif

/*
 * Level 2 BLAS
 */

#if MANUFACTURE != CRAY
  /*
   * Single precision
   */

  #define            PLA_sgemv(           trans,       m, n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    sgemv  ) (       trans,       m, n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_ssymv(     uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    ssymv  ) ( uplo,                 n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_strmv(     uplo, trans, diag,    n,        a, lda, x, incx ) \
          PLA_FOR2C(     strmv  ) ( uplo, trans, diag,    n,        a, lda, x, incx, 1, 1, 1 )
  #define            PLA_strsv(     uplo, trans, diag,    n,        a, lda, x, incx ) \
          PLA_FOR2C(     strsv  ) ( uplo, trans, diag,    n,        a, lda, x, incx, 1, 1, 1 )

  #define           PLA_sger(            m, n, alpha, x, incx, y, incy, a, lda ) \
          PLA_FOR2C(    sger ) (         m, n, alpha, x, incx, y, incy, a, lda )
  #define           PLA_ssyr(      uplo,    n, alpha, x, incx,          a, lda ) \
          PLA_FOR2C(    ssyr   ) ( uplo,    n, alpha, x, incx,          a, lda, 1 )
  #define           PLA_ssyr2(     uplo,    n, alpha, x, incx, y, incy, a, lda ) \
          PLA_FOR2C(    ssyr2  ) ( uplo,    n, alpha, x, incx, y, incy, a, lda, 1 )


  /*
   * Double precision
   */

  #define            PLA_dgemv(           trans,       m, n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    dgemv  ) (       trans,       m, n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_dsymv(     uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    dsymv  ) ( uplo,                 n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_dtrmv(     uplo, trans, diag,    n,        a, lda, x, incx ) \
          PLA_FOR2C(     dtrmv  ) ( uplo, trans, diag,    n,        a, lda, x, incx, 1, 1, 1 )
  #define            PLA_dtrsv(     uplo, trans, diag,    n,        a, lda, x, incx ) \
          PLA_FOR2C(     dtrsv  ) ( uplo, trans, diag,    n,        a, lda, x, incx, 1, 1, 1 )

  #define           PLA_dger(            m, n, alpha, x, incx, y, incy, a, lda ) \
          PLA_FOR2C(    dger ) (         m, n, alpha, x, incx, y, incy, a, lda )
  #define           PLA_dsyr(      uplo,    n, alpha, x, incx,          a, lda ) \
          PLA_FOR2C(    dsyr   ) ( uplo,    n, alpha, x, incx,          a, lda, 1 )
  #define           PLA_dsyr2(     uplo,    n, alpha, x, incx, y, incy, a, lda ) \
          PLA_FOR2C(    dsyr2  ) ( uplo,    n, alpha, x, incx, y, incy, a, lda, 1 )

  /*
   * Single precision complex
   */

  #define            PLA_cgemv(           trans,       m, n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    cgemv  ) (       trans,       m, n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_csymv(     uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    csymv  ) ( uplo,                 n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_chemv(     uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    chemv  ) ( uplo,                 n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_ctrmv(     uplo, trans, diag,    n,        a, lda, x, incx ) \
          PLA_FOR2C(     ctrmv  ) ( uplo, trans, diag,    n,        a, lda, x, incx, 1, 1, 1 )
  #define            PLA_ctrsv(     uplo, trans, diag,    n,        a, lda, x, incx ) \
          PLA_FOR2C(     ctrsv  ) ( uplo, trans, diag,    n,        a, lda, x, incx, 1, 1, 1 )

  #define           PLA_cger(            m, n, alpha, x, incx, y, incy, a, lda ) \
          PLA_FOR2C(    cgeru ) (         m, n, alpha, x, incx, y, incy, a, lda )
  #define           PLA_cher(      uplo,    n, alpha, x, incx,          a, lda ) \
          PLA_FOR2C(    cher   ) ( uplo,    n, alpha, x, incx,          a, lda, 1 )
  #define           PLA_cher2(     uplo,    n, alpha, x, incx, y, incy, a, lda ) \
          PLA_FOR2C(    cher2  ) ( uplo,    n, alpha, x, incx, y, incy, a, lda, 1 )

  /*
   * Double precision complex
   */
  #define            PLA_zgemv(           trans,       m, n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    zgemv  ) (       trans,       m, n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_zsymv(     uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    zsymv  ) ( uplo,                 n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_zhemv(     uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
          PLA_FOR2C (    zhemv  ) ( uplo,                 n, alpha, a, lda, x, incx, beta, y, incy, 1 )
  #define            PLA_ztrmv(     uplo, trans, diag,    n,        a, lda, x, incx ) \
          PLA_FOR2C(     ztrmv  ) ( uplo, trans, diag,    n,        a, lda, x, incx, 1, 1, 1 )
  #define            PLA_ztrsv(     uplo, trans, diag,    n,        a, lda, x, incx ) \
          PLA_FOR2C(     ztrsv  ) ( uplo, trans, diag,    n,        a, lda, x, incx, 1, 1, 1 )

  #define           PLA_zger(            m, n, alpha, x, incx, y, incy, a, lda ) \
          PLA_FOR2C(    zgeru ) (         m, n, alpha, x, incx, y, incy, a, lda )
  #define           PLA_zher(      uplo,    n, alpha, x, incx,          a, lda ) \
          PLA_FOR2C(    zher   ) ( uplo,    n, alpha, x, incx,          a, lda, 1 )
  #define           PLA_zher2(     uplo,    n, alpha, x, incx, y, incy, a, lda ) \
          PLA_FOR2C(    zher2  ) ( uplo,    n, alpha, x, incx, y, incy, a, lda, 1 )
 

#else /*  MANUFACTURE == CRAY */

  /*
   * Single precision
   */

  #define PLA_sgemv(                              trans,       m, n, alpha, a, lda, x, incx, beta, y, incy ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_ssymv(                        uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_strmv(                        uplo, trans, diag,    n,        a, lda, x, incx ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_strsv(                        uplo, trans, diag,    n,        a, lda, x, incx ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_sger(                            m, n, alpha, x, incx, y, incy, a, lda ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_ssyr(                      uplo,    n, alpha, x, incx,          a, lda ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_ssyr2(                     uplo,    n, alpha, x, incx, y, incy, a, lda ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  /*
   * Double precision
   */

  #define PLA_dgemv(                              trans,       m, n, alpha, a, lda, x, incx, beta, y, incy ) \
             SGEMV  ( PLA_CP2FCD ( trans ),                                                                  \
                                                              m, n, alpha, a, lda, x, incx, beta, y, incy )
  #define PLA_dsymv(                        uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
             SSYMV  ( PLA_CP2FCD ( uplo  ),                                                                  \
                                                                 n, alpha, a, lda, x, incx, beta, y, incy )

  #define PLA_dtrmv(                        uplo, trans, diag,    n,        a, lda, x, incx ) \
             STRMV  ( PLA_CP2FCD ( uplo  ),                                                   \
                      PLA_CP2FCD ( trans ),                                                   \
                      PLA_CP2FCD ( diag  ),                                                   \
                                                                 n,        a, lda, x, incx )
  #define PLA_dtrsv(                        uplo, trans, diag,    n,        a, lda, x, incx ) \
             STRSV  ( PLA_CP2FCD ( uplo  ),                                                   \
                      PLA_CP2FCD ( trans ),                                                   \
                      PLA_CP2FCD ( diag  ),                                                   \
                                                                 n,        a, lda, x, incx )

  #define PLA_dger(                            m, n, alpha, x, incx, y, incy, a, lda ) \
             SGER (                           m, n, alpha, x, incx, y, incy, a, lda )
  #define PLA_dsyr(                      uplo,    n, alpha, x, incx,          a, lda ) \
             SSYR ( PLA_CP2FCD ( uplo ),                                               \
                                                 n, alpha, x, incx,          a, lda )
  #define PLA_dsyr2(                     uplo,    n, alpha, x, incx, y, incy, a, lda ) \
             SSYR2 ( PLA_CP2FCD ( uplo ), n, alpha, x, incx, y, incy, a, lda )

  /* prototypes */

  void  SGEMV( _fcd , int *, int *, double *, double *, int *, double *, int *, 
	       double *, double *, int * );
  void  SSYMV( _fcd , int *, double *, double *, int *, double *, int *, 
	       double *, double *, int * );
  void  STRMV( _fcd , _fcd, _fcd, int *, double *, int *, double *, int * ); 
  void  STRSV( _fcd , _fcd, _fcd, int *, double *, int *, double *, int * ); 

  void  SGER( int *, int *, double *, double *, int *, double *, int *, double *, int * ); 
  void  SSYR( _fcd,  int *, double *, double *, int *, double *, int * ); 
  void  SSYR2(_fcd,  int *, double *, double *, int *, double *, int *, double *, int * ); 



  /*
   * Single precision complex
   */

  #define PLA_cgemv(                              trans,       m, n, alpha, a, lda, x, incx, beta, y, incy ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_csymv(                        uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_chemv(                        uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_ctrmv(                        uplo, trans, diag,    n,        a, lda, x, incx ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_ctrsv(                        uplo, trans, diag,    n,        a, lda, x, incx ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_cger(                            m, n, alpha, x, incx, y, incy, a, lda ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_csyr(                      uplo,    n, alpha, x, incx,          a, lda ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )
  #define PLA_csyr2(                     uplo,    n, alpha, x, incx, y, incy, a, lda ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  /*
   * Double precision complex
   */

  #define PLA_zgemv(                              trans,       m, n, alpha, a, lda, x, incx, beta, y, incy ) \
             CGEMV  ( PLA_CP2FCD ( trans ),                                                                  \
                                                              m, n, alpha, a, lda, x, incx, beta, y, incy )
  #define PLA_zsymv(                        uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
             CSYMV  ( PLA_CP2FCD ( uplo  ),                                                                  \
                                                                 n, alpha, a, lda, x, incx, beta, y, incy )
  #define PLA_zhemv(                        uplo,                 n, alpha, a, lda, x, incx, beta, y, incy ) \
             CHEMV  ( PLA_CP2FCD ( uplo  ),                                                                  \
                                                                 n, alpha, a, lda, x, incx, beta, y, incy )
  #define PLA_ztrmv(                        uplo, trans, diag,    n,        a, lda, x, incx ) \
             CTRMV  ( PLA_CP2FCD ( uplo  ),                                                   \
                      PLA_CP2FCD ( trans ),                                                   \
                      PLA_CP2FCD ( diag  ),                                                   \
                                                                 n,        a, lda, x, incx )
  #define PLA_ztrsv(                        uplo, trans, diag,    n,        a, lda, x, incx ) \
             CTRSV  ( PLA_CP2FCD ( uplo  ),                                                   \
                      PLA_CP2FCD ( trans ),                                                   \
                      PLA_CP2FCD ( diag  ),                                                   \
                                                                 n,        a, lda, x, incx )

  #define PLA_zger(                            m, n, alpha, x, incx, y, incy, a, lda ) \
             CGERU (                           m, n, alpha, x, incx, y, incy, a, lda )
  #define PLA_zsyr(                      uplo,    n, alpha, x, incx,          a, lda ) \
             CSYR ( PLA_CP2FCD ( uplo ),                                               \
                                                 n, alpha, x, incx,          a, lda )
  #define PLA_zsyr2(                     uplo,    n, alpha, x, incx, y, incy, a, lda ) \
             CSYR2 ( PLA_CP2FCD ( uplo ),                                              \
                                                 n, alpha, x, incx, y, incy, a, lda )

  /* prototypes */

  void  CGEMV( _fcd , int *, int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, 
	       PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  void  CSYMV( _fcd , int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, 
	       PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  void  CTRMV( _fcd , _fcd, _fcd, int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * ); 
  void  CTRSV( _fcd , _fcd, _fcd, int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * ); 

  void  CGER( int *, int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * ); 
  void  CSYR( _fcd,  int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * ); 
  void  CSYR2( _fcd,  int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * ); 
#endif

/*
 * Level 3 BLAS
 */

#if MANUFACTURE != CRAY

  /*
   * Single precision
   */

  #define           PLA_sgemm( transa, transb,  m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    sgemm  ) (   transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_ssymm(     side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    ssymm  ) ( side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_ssyrk(           uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc )  \
          PLA_FOR2C(    ssyrk  ) (       uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc, 1, 1 )
  #define           PLA_ssyr2k(          uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    ssyr2k ) (       uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_strmm(     side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               )  \
          PLA_FOR2C(    strmm  ) ( side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb, 1, 1, 1, 1               )
  #define           PLA_strsm(     side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               )  \
          PLA_FOR2C(    strsm  ) ( side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb, 1, 1, 1, 1               )

  /*
   * Double precision
   */

  #define           PLA_dgemm(                 transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    dgemm  ) (             transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_dsymm(     side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    dsymm  ) ( side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_dsyrk(           uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc )  \
          PLA_FOR2C(    dsyrk  ) (       uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc, 1, 1 )
  #define           PLA_dsyr2k(          uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    dsyr2k ) (       uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_dtrmm(     side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               )  \
          PLA_FOR2C(    dtrmm  ) ( side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb, 1, 1, 1, 1               )
  #define           PLA_dtrsm(     side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               )  \
          PLA_FOR2C(    dtrsm  ) ( side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb, 1, 1, 1, 1               )

  /*
   * Single precision complex
   */

  #define           PLA_cgemm(                 transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    cgemm  ) (             transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_csymm(     side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    csymm  ) ( side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_chemm(     side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    chemm  ) ( side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_csyrk(           uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc )  \
          PLA_FOR2C(    csyrk  ) (       uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc, 1, 1 )
  #define           PLA_cherk(           uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc )  \
          PLA_FOR2C(    cherk  ) (       uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc, 1, 1 )
  #define           PLA_csyr2k(          uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    csyr2k ) (       uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_cher2k(          uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    cher2k ) (       uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_ctrmm(     side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               )  \
          PLA_FOR2C(    ctrmm  ) ( side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb, 1, 1, 1, 1               )
  #define           PLA_ctrsm(     side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               )  \
          PLA_FOR2C(    ctrsm  ) ( side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb, 1, 1, 1, 1               )

  /*
   * Double precision complex
   */
  #define           PLA_zgemm(                 transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    zgemm  ) (             transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_zsymm(     side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    zsymm  ) ( side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_zhemm(     side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    zhemm  ) ( side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_zsyrk(           uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc )  \
          PLA_FOR2C(    zsyrk  ) (       uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc, 1, 1 )
  #define           PLA_zherk(           uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc )  \
          PLA_FOR2C(    zherk  ) (       uplo, trans,                   n, k, alpha, a, lda,         beta, c, ldc, 1, 1 )
  #define           PLA_zsyr2k(          uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    zsyr2k ) (       uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_zher2k(          uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc )  \
          PLA_FOR2C(    zher2k ) (       uplo, trans,                   n, k, alpha, a, lda, b, ldb, beta, c, ldc, 1, 1 )
  #define           PLA_ztrmm(     side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               )  \
          PLA_FOR2C(    ztrmm  ) ( side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb, 1, 1, 1, 1               )
  #define           PLA_ztrsm(     side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               )  \
          PLA_FOR2C(    ztrsm  ) ( side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb, 1, 1, 1, 1               )

#else /* MANUFACTURE == CRAY */


  /*
   * Single precision
   */


  #define PLA_sgemm(                           transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_ssymm(               side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_ssyrk(                    uplo, trans,                    n, k, alpha, a, lda,         beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_ssyr2k(                   uplo, trans,                    n, k, alpha, a, lda, b, ldb, beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_strmm(               side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_strsm(               side, uplo, transa,         diag, m, n,     alpha, a, lda, b, ldb              ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  /*
   * Double precision
   */
  #define PLA_dgemm( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc ) \
             SGEMM  ( PLA_CP2FCD ( transa ),   \
                      PLA_CP2FCD ( transb ),   \
                      m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
  #define PLA_dsymm(               side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc ) \
             SSYMM  ( PLA_CP2FCD ( side   ),                                                                        \
                      PLA_CP2FCD ( uplo   ),                                                                        \
                                                                    m, n,    alpha, a, lda, b, ldb, beta, c, ldc )
  #define PLA_dsyrk(                    uplo, trans,                    n, k, alpha, a, lda,         beta, c, ldc ) \
             SSYRK  ( PLA_CP2FCD ( uplo   ),                                                                        \
                      PLA_CP2FCD ( trans  ),                                                                        \
                                                                       n, k, alpha, a, lda,         beta, c, ldc )
  #define PLA_dsyr2k(                   uplo, trans,                    n, k, alpha, a, lda, b, ldb, beta, c, ldc ) \
             SSYR2K ( PLA_CP2FCD ( uplo   ),                                                                        \
                      PLA_CP2FCD ( trans  ),                                                                        \
                                                                       n, k, alpha, a, lda, b, ldb, beta, c, ldc )
  #define PLA_dtrmm(               side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               ) \
             STRMM  ( PLA_CP2FCD ( side   ),                                                                        \
                      PLA_CP2FCD ( uplo   ),                                                                        \
                      PLA_CP2FCD ( transa ),                                                                        \
                      PLA_CP2FCD ( diag   ),                                                                        \
                                                                   m, n,     alpha, a, lda, b, ldb               )
  #define PLA_dtrsm(               side, uplo, transa,        diag, m, n,     alpha, a, lda, b, ldb               ) \
             STRSM  ( PLA_CP2FCD ( side   ),                                                                        \
                      PLA_CP2FCD ( uplo   ),                                                                        \
                      PLA_CP2FCD ( transa ),                                                                        \
                      PLA_CP2FCD ( diag   ),                                                                        \
                                                                   m, n,     alpha, a, lda, b, ldb               )

  /* prototypes */

  void  SGEMM( _fcd , _fcd, int *, int *, int *, double *, double *, int *, double *, int *, 
	       double *, double *, int * );
  void  SSYMM( _fcd , _fcd, int *, int *,        double *, double *, int *, double *, int *, 
	       double *, double *, int * );
  void  SSYRK( _fcd , _fcd, int *, int *,        double *, double *, int *, 
	       double *, double *, int * );
  void  SSYR2K(_fcd , _fcd, int *, int *,        double *, double *, int *, double *, int *,
	       double *, double *, int * );
  void  STRMM( _fcd , _fcd, _fcd, _fcd, int *, int *, double *, double *, int *, double *, int * );
  void  STRSM( _fcd , _fcd, _fcd, _fcd, int *, int *, double *, double *, int *, double *, int * );


   /*
    * Single precision complex
    */

  #define PLA_cgemm(                           transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_csymm(               side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_cher2(               uplo,                          n, alpha, x, incx, y,incy,  a, lda ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_cher(               uplo,                          n, alpha, x, incx,          a, lda ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_csyrk(                    uplo, trans,                    n, k, alpha, a, lda,         beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_cherk(                    uplo, trans,                    n, k, alpha, a, lda,         beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_csyr2k(                   uplo, trans,                    n, k, alpha, a, lda, b, ldb, beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_cher2k(                   uplo, trans,                    n, k, alpha, a, lda, b, ldb, beta, c, ldc ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_ctrmm(               side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  #define PLA_ctrsm(               side, uplo, transa,        diag, m, n,     alpha, a, lda, b, ldb               ) \
                       PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

  /*
   * double precision complex
   */
  #define PLA_zgemm(                           transa, transb,       m, n, k, alpha, a, lda, b, ldb, beta, c, ldc ) \
             CGEMM  ( PLA_CP2FCD ( transa ),                                                                        \
                      PLA_CP2FCD ( transb ),                                                                        \
                                                                    m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
  #define PLA_zsymm(               side, uplo,                       m, n,    alpha, a, lda, b, ldb, beta, c, ldc ) \
             CSYMM  ( PLA_CP2FCD ( side   ),                                                                        \
                      PLA_CP2FCD ( uplo   ),                                                                        \
                                                                    m, n,    alpha, a, lda, b, ldb, beta, c, ldc )
  #define PLA_zher(               uplo,                        n, alpha, x, incx,          a, lda ) \
             CHER   (                                                                       \
                      PLA_CP2FCD ( uplo   ),                                                                        \
                                                                     n, alpha, x, incx,          a, lda )
  #define PLA_zsyrk(                    uplo, trans,                    n, k, alpha, a, lda,         beta, c, ldc ) \
             CSYRK  ( PLA_CP2FCD ( uplo   ),                                                                        \
                      PLA_CP2FCD ( trans  ),                                                                        \
                                                                       n, k, alpha, a, lda,         beta, c, ldc )
  #define PLA_zherk(                    uplo, trans,                    n, k, alpha, a, lda,         beta, c, ldc ) \
             CHERK  ( PLA_CP2FCD ( uplo   ),                                                                        \
                      PLA_CP2FCD ( trans  ),                                                                        \
                                                                       n, k, alpha, a, lda,         beta, c, ldc )
  #define PLA_zsyr2k(                   uplo, trans,                    n, k, alpha, a, lda, b, ldb, beta, c, ldc ) \
             CSYR2K ( PLA_CP2FCD ( uplo   ),                                                                        \
                      PLA_CP2FCD ( trans  ),                                                                        \
                                                                       n, k, alpha, a, lda, b, ldb, beta, c, ldc )
  #define PLA_zher2(                   uplo, n, alpha, x, incx, y, incy, a, lda        ) \
             CHER2K ( PLA_CP2FCD ( uplo   ),                                                                        \
                                                                       n, alpha, x, incx, y, incy, a, lda        )
  #define PLA_ztrmm(               side, uplo, transa,         diag, m, n,    alpha, a, lda, b, ldb               ) \
             CTRMM  ( PLA_CP2FCD ( side   ),                                                                        \
                      PLA_CP2FCD ( uplo   ),                                                                        \
                      PLA_CP2FCD ( transa ),                                                                        \
                      PLA_CP2FCD ( diag   ),                                                                        \
                                                                   m, n,     alpha, a, lda, b, ldb               )
  #define PLA_ztrsm(               side, uplo, transa,        diag, m, n,     alpha, a, lda, b, ldb               ) \
             CTRSM  ( PLA_CP2FCD ( side   ),                                                                        \
                      PLA_CP2FCD ( uplo   ),                                                                        \
                      PLA_CP2FCD ( transa ),                                                                        \
                      PLA_CP2FCD ( diag   ),                                                                        \
                                                                   m, n,     alpha, a, lda, b, ldb               )

  /* prototypes */

  void  CGEMM( _fcd , _fcd, int *, int *, int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, 
	       PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  void  CSYMM( _fcd , _fcd, int *, int *,        PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, 
	       PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  void  CHEMM( _fcd , _fcd, int *, int *,        PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, 
	       PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  void  CSYRK( _fcd , _fcd, int *, int *,        PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, 
	       PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  void  CHERK( _fcd , _fcd, int *, int *,        PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, 
	       PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  void  CSYR2K(_fcd , _fcd, int *, int *,        PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *,
	       PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  void  CHER2K(_fcd , _fcd, int *, int *,        PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *,
	       PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
  void  CTRMM( _fcd , _fcd, _fcd, _fcd, int *, int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
  void  CTRSM( _fcd , _fcd, _fcd, _fcd, int *, int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );

#endif

/*
 * LAPACK interface
 */

#if MANUFACTURE != CRAY
#define PLA_dsteqr( compz, n, d, e, z, ldz, work, info ) \
          PLA_FOR2C( dsteqr  ) ( compz, n, d, e, z, ldz, work, info, 1 )
#else
#define PLA_dsteqr( compz, n, d, e, z, ldz, work, info ) \
          SSTEQR  ( PLA_CP2FCD ( compz ), n, d, e, z, ldz, work, info )
#endif

#if MANUFACTURE != CRAY
#if MANUFACTURE != PC
#define PLA_dsteqr_x( compz, m, n, d, e, z, ldz, work, info ) \
          PLA_FOR2C( dsteqr_x )( compz, m, n, d, e, z, ldz, work, info, 1 )

/* LAPACK Cholesky factorization */
#define PLA_spotrf( uplo, n, a, lda, info ) \
          PLA_FOR2C(    spotrf  ) (   uplo, n, a, lda, info, 1 )

#define PLA_dpotrf( uplo, n, a, lda, info ) \
          PLA_FOR2C(    dpotrf  ) (   uplo, n, a, lda, info, 1 )

#define PLA_cpotrf( uplo, n, a, lda, info ) \
          PLA_FOR2C(    cpotrf  ) (   uplo, n, a, lda, info, 1 )

#define PLA_zpotrf( uplo, n, a, lda, info ) \
          PLA_FOR2C(    zpotrf  ) (   uplo, n, a, lda, info, 1 )
#else
#define PLA_dsteqr_x( compz, m, n, d, e, z, ldz, work, info ) \
          PLA_FOR2C( dsteqr_x )( compz, m, n, d, e, z, ldz, work, info, 1 )

#define PLA_dbdsqr_x( uplo, n, nv, nru, ncc, d, e, v, ldv, \
		      u, ldu, c, ldc, work, info ) \
          PLA_FOR2C( dbdsqr_x)( uplo, n, nv, nru, ncc, d, e, v, ldv, \
		      u, ldu, c, ldc, work, info, 1 ) 

/* LAPACK Cholesky factorization */
#define PLA_spotrf( uplo, n, a, lda, info ) \
          PLA_FOR2C( spotrf_ )(   uplo, n, a, lda, info, 1 )

#define PLA_dpotrf( uplo, n, a, lda, info ) \
          PLA_FOR2C( dpotrf_ )(   uplo, n, a, lda, info, 1 )

#define PLA_cpotrf( uplo, n, a, lda, info ) \
          PLA_FOR2C( cpotrf_ )(   uplo, n, a, lda, info, 1 )

#define PLA_zpotrf( uplo, n, a, lda, info ) \
          PLA_FOR2C( zpotrf_ )(   uplo, n, a, lda, info, 1 )
#endif
#else
#define PLA_dsteqr_x( compz, m, n, d, e, z, ldz, work, info ) \
          SSTEQR_X  ( PLA_CP2FCD ( compz ), m, n, d, e, z, ldz, work, info )

/* LAPACK Cholesky factorization */
#define PLA_spotrf( uplo, n, a, lda, info ) \
            PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

#define PLA_dpotrf( uplo, n, a, lda, info ) \
             SPOTRF ( PLA_CP2FCD ( uplo ), n, a, lda, info )

#define PLA_cpotrf( uplo, n, a, lda, info ) \
            PLA_Abort( "32 bit fl. point not supported ", __LINE__, __FILE__ )

#define PLA_zpotrf( uplo, n, a, lda, info ) \
             CPOTRF ( PLA_CP2FCD ( uplo ), n, a, lda, info )

#endif

