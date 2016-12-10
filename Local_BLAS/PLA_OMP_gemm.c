/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/


#include <stdlib.h>
#include <PLA.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int getenv_int(const char *varname) {
    char *p_env = getenv(varname);
    if (p_env == NULL) return 0;
    return atoi(p_env);
}

/*----------------------------------------------------------------------
 *  PLA_OMP_SGEMM - Single Precision General Matrix Matrix Multiply
 *--------------------------------------------------------------------*/
void PLA_OMP_sgemm (
          char*  transa, 
	  char*  transb,
	  int m, int n, int k, 
	  float* alpha, 
	  float* a, int lda, 
	  float* b, int ldb, 
	  float* beta, 
	  float* c, int ldc )
{
    int nt, rows, cols, tid;

    /* Executable Statements */

    /* printf("Entry to PLA_OMP_sgemm\n"); */

    /* Determine how many threads we can get */

#ifdef _OPENMP
    nt = omp_get_max_threads();
    /* printf("omp_get_max_threads returned nt = %d\n", nt ); */
#else
    if ((nt = getenv_int("OMP_NUM_THREADS")) == 0) nt = 1;
    /* printf("omp disabled.  num threads statically set to %d\n", nt ); */
#endif

    /* Make sure there's enough computations to warrant multithreading ... */

    if ( nt <= 1 || n < nt || n*m < 100 ) { 
	PLA_sgemm(
	        transa, transb,
		&m, &n, &k,
		alpha, 
		a, &lda, 
		b, &ldb, 
	 	beta, 
		c, &ldc); 
    } else {

        /* Use OpenMP to parallelize SGEMM computation */

        if ((rows = getenv_int("PLA_OMP_ROWS")) == 0) rows = 1;
	if ((cols = getenv_int("PLA_OMP_COLS")) == 0) cols = nt;

	if ( rows*cols != nt ) {
	  printf("PLA_OMP_gemm: Error: Rows and cols don't match!!\n");
	  exit(1);
	}

        /*
	 * We decompose C, A, and B into blocks of rows and columns as
	 * specified by the rows and cols determined above.
	 */

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static,1)
#endif
        for (tid=0; tid < nt; tid++) {

	    int zbrow, zbcol;     /* zero based row & col for this thread */
	    int mb, nb;           /* number of rows and cols in a typical 
				     block (not in last row or col)  */
	    int ib, jb;           /* number of rows and cols in the block 
				     for this thread */

	    float * ablock, * bblock; /* blocks of A and B to use */

	    /* Use the thread id to determine the block of the matrix C
	       that will be computed by this thread.  Threads are assigned 
	       blocks of C in column major order.  Thus ... */

	    zbrow = tid % rows;
	    zbcol = tid / rows;

	    /* The upper left corner (of C) of the block for this thread
	       is C( zbrow*mb, zbcol*nb ).  Now we compute the size of the
	       typical block. */

	    mb = (m+rows-1) / rows;  /* mb = ceiling( m / rows ) */
	    nb = (n+cols-1) / cols;  /* nb = ceiling( n / cols ) */

	    /* If the matrix dimensions don't divide evenly into our
	       rows and columns, we make up for it in the last row/col */

            ib = min( mb, m-(zbrow*mb) );  
            jb = min( nb, n-(zbcol*nb) );

	    /* printf("Tid=%d, ib=%d, jb=%d \n", tid, ib, jb); */

	    /* Now get ptrs to blocks of A and B to pass in to the BLAS */

	    if ( *transa == 'N' )
		ablock = &a[ (zbrow*mb) ];     /* element A[ zbrow*mb +1 , 1 ] */
	    else
		ablock = &a[ (zbrow*mb) * lda ];  /* element A[ 1, zbrow*mb +1 ] */

	    if ( *transb == 'N' )
		bblock = &b[ (zbcol*nb) * ldb ];  /* element B[ 1, zbcol*nb + 1 ] */
	    else
		bblock = &b[ (zbcol*nb) ];       /* element B[ zbcol*nb + 1, 1 ] */

	    PLA_sgemm(
		transa, transb,
		&ib, &jb, &k,
		alpha, 
		ablock, &lda, 
		bblock, &ldb, 
		beta, 
		&c[ (zbcol*nb) * ldc + (zbrow*mb)], &ldc);  /* element C[ zbrow*mb+1, zbcol*nb+1 ] */
	}
        /* End of parallel section */
    }

    /* printf("Return from PLA_OMP_sgemm\n"); */

    /* End of PLA_OMP_SGEMM */
}

/*----------------------------------------------------------------------
 *  PLA_OMP_DGEMM - Double Precision General Matrix Matrix Multiply
 *--------------------------------------------------------------------*/
void PLA_OMP_dgemm (
          char*  transa, 
	  char*  transb,
	  int m, int n, int k, 
	  double* alpha, 
	  double* a, int lda, 
	  double* b, int ldb, 
	  double* beta, 
	  double* c, int ldc )
{
    int nt, rows, cols, tid;

    /* Executable Statements */

    /* printf("Entry to PLA_OMP_dgemm\n"); */

    /* Determine how many threads we can get */

#ifdef _OPENMP
    nt = omp_get_max_threads();
    /* printf("omp_get_max_threads returned nt = %d\n", nt ); */
#else
    if ((nt = getenv_int("OMP_NUM_THREADS")) == 0) nt = 1;
    /* printf("omp disabled.  num threads statically set to %d\n", nt ); */
#endif

    /* Make sure there's enough computations to warrant multithreading ... */

    if ( nt <= 1 || n < nt || n*m < 1 ) { 
	PLA_dgemm(
	        transa, transb,
		&m, &n, &k,
		alpha, 
		a, &lda, 
		b, &ldb, 
	 	beta, 
		c, &ldc); 
    } else {

        /* Use OpenMP to parallelize DGEMM computation */

        if ((rows = getenv_int("PLA_OMP_ROWS")) == 0) rows = 1;
	if ((cols = getenv_int("PLA_OMP_COLS")) == 0) cols = nt;

	if ( rows*cols != nt ) {
	  printf("PLA_OMP_gemm: Error: Rows and cols don't match!!\n");
	  exit(1);
	}

        /*
	 * We decompose C, A, and B into blocks of rows and columns as
	 * specified by the rows and cols determined above.
	 */

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static,1)
#endif
        for (tid=0; tid < nt; tid++) {

	    int zbrow, zbcol;     /* zero based row & col for this thread */
	    int mb, nb;           /* number of rows and cols in a typical 
				     block (not in last row or col)  */
	    int ib, jb;           /* number of rows and cols in the block 
				     for this thread */

	    double * ablock, * bblock; /* blocks of A and B to use */

	    /* Use the thread id to determine the block of the matrix C
	       that will be computed by this thread.  Threads are assigned 
	       blocks of C in column major order.  Thus ... */

	    zbrow = tid % rows;
	    zbcol = tid / rows;

	    /* The upper left corner (of C) of the block for this thread
	       is C( zbrow*mb, zbcol*nb ).  Now we compute the size of the
	       typical block. */

	    mb = (m+rows-1) / rows;  /* mb = ceiling( m / rows ) */
	    nb = (n+cols-1) / cols;  /* nb = ceiling( n / cols ) */

	    /* If the matrix dimensions don't divide evenly into our
	       rows and columns, we make up for it in the last row/col */

            ib = min( mb, m-(zbrow*mb) );  
            jb = min( nb, n-(zbcol*nb) );

	    /* printf("Tid=%d, ib=%d, jb=%d \n", tid, ib, jb); */

	    /* Now get ptrs to blocks of A and B to pass in to the BLAS */

	    if ( *transa == 'N' )
		ablock = &a[ (zbrow*mb) ];     /* element A[ zbrow*mb +1 , 1 ] */
	    else
		ablock = &a[ (zbrow*mb) * lda ];  /* element A[ 1, zbrow*mb +1 ] */

	    if ( *transb == 'N' )
		bblock = &b[ (zbcol*nb) * ldb ];  /* element B[ 1, zbcol*nb + 1 ] */
	    else
		bblock = &b[ (zbcol*nb) ];       /* element B[ zbcol*nb + 1, 1 ] */

	    PLA_dgemm(
		transa, transb,
		&ib, &jb, &k,
		alpha, 
		ablock, &lda, 
		bblock, &ldb, 
		beta, 
		&c[ (zbcol*nb) * ldc + (zbrow*mb)], &ldc);  /* element C[ zbrow*mb+1, zbcol*nb+1 ] */
	}
        /* End of parallel section */
    }

    /* printf("Return from PLA_OMP_dgemm\n"); */

}
