C
C   PLAPACK General Linear Solver example application.
C
C   Purpose:   Implement the solution of a problem requiring
C              a general linear solver, including the problem
C              generation phase, the solution of the system,
C              and the retrieval of the solution data.
C
C   Datatypes: MPI_FLOAT, MPI_DOUBLE, MPI_DOUBLE_COMPLEX 
C               
C


      PROGRAM MAIN

      implicit none

      include "mpif.h"
      include "PLAf.h"
C  
C MPI_Comm
C     
      INTEGER COMM
C
C PLA_Template
C    
      INTEGER TEMPL
C
C PLA_Obj
C
      INTEGER 
     &     A,
     &     PIVOTS,
     &     X,
     &     B,
     &     B_NORM,
     &     INDEX,
     &     MINUS_ONE,
     &     ZERO,
     &     ONE
C
      DOUBLE PRECISION 
     &     OPERATION_COUNT, 
     &     B_NORM_VALUE, 
     &     TIME
      INTEGER 
     &     IERROR,
     &     SIZE, 
     &     NB_DISTR,
     &     ME,
     &     NPROCS,
     &     NPROWS,
     &     NPCOLS,
     &     DUMMY,
     &     INFO 
C
C MPI_Datatype 
C
      INTEGER DATATYPE


      COMM = MPI_COMM_NULL
      TEMPL = 0
      A = 0 
      PIVOTS = 0 
      X = 0 
      B = 0 
      B_NORM = 0 
      INDEX = 0
      MINUS_ONE = 0
      ZERO = 0
      ONE = 0

      CALL MPI_INIT( IERROR )

C
C Get problem size and distribution block size and broadcast 
C
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, ME, IERROR )
      if ( ME .EQ. 0 ) THEN
         PRINT *, "-- enter processor mesh dimension ( rows cols ):"
         READ  *, NPROWS, NPCOLS
         PRINT *, "-- enter matrix size, distr. block size:"
         READ  *, SIZE, NB_DISTR
      ENDIF

      CALL MPI_BCAST( NPROWS,   1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &     IERROR )
      CALL MPI_BCAST( NPCOLS,   1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &     IERROR )
      CALL MPI_BCAST( SIZE,     1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &     IERROR )
      CALL MPI_BCAST( NB_DISTR, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     &     IERROR )

C
C Two different ways to create a 2D communicator 
C
C  PLA_Comm1dTo2dRatio_F ( MPI_COMM_WORLD, 0.25d00, COMM )
      CALL PLA_Comm1dTo2d_F( MPI_COMM_WORLD, NPROWS, NPCOLS, 
     &     COMM, IERROR )
C
C Initialize PLAPACK 
C
      CALL PLA_INIT_F( COMM, IERROR )
C
C Create object distribution template 
C
      CALL PLA_TempCreate_F( NB_DISTR, 0, TEMPL, IERROR )
C
C Set the datatype : MPI_DOUBLE, MPI_FLOAT or MPI_DOUBLE_COMPLEX 
C
      DATATYPE = MPI_DOUBLE_PRECISION
C
C Create objects for problem to be solved 
C
      CALL PLA_MatrixCreate_F(  DATATYPE, SIZE, SIZE, TEMPL,
     &                           PLA_ALIGN_FIRST, PLA_ALIGN_FIRST, 
     &                           A, IERROR )

      CALL PLA_MvectorCreate_F( DATATYPE, SIZE, 1, 
     &                           TEMPL, PLA_ALIGN_FIRST, 
     &                           X, IERROR )

      CALL PLA_MvectorCreate_F( DATATYPE, SIZE, 1, 
     &                           TEMPL, PLA_ALIGN_FIRST, 
     &                           B, IERROR )

      CALL PLA_MvectorCreate_F( MPI_INTEGER, 
     &                           SIZE, 1, TEMPL, PLA_ALIGN_FIRST,
     &                           PIVOTS, IERROR )
C
C Create duplicated scalar constants with same datatype and template as A 
C
      CALL PLA_CreateConstantsConfTo_F( A, MINUS_ONE, ZERO, ONE, 
     &     IERROR )
C
C Create a problem to be solved: A x = b 
C
      CALL CREATE_PROBLEM_F( A, X, B )
C
C Start timing 
      CALL MPI_BARRIER( MPI_COMM_WORLD, IERROR )
      TIME = MPI_WTIME( )
C
C Factor P A -> L U overwriting lower triangular portion of A with L, upper, U 
C
      CALL PLA_LU_F( A, PIVOTS, INFO )

      IF ( INFO .NE. 0 ) THEN
         PRINT *, "Zero pivot encountered at row", INFO
      ELSE
C
C  Apply the permutations to the right hand sides 
C
         CALL PLA_ApplyPivotsToRows_F ( B, PIVOTS, IERROR )
C
C Solve L y = b, overwriting b with y 
C
         CALL PLA_Trsv_F( PLA_LOWER_TRIANGULAR, PLA_NO_TRANSPOSE, 
     &                    PLA_UNIT_DIAG, A, B, IERROR )
C
C Solve U x = y (=b), overwriting b with x
C
         CALL PLA_Trsv_F( PLA_UPPER_TRIANGULAR, PLA_NO_TRANSPOSE,
     &                    PLA_NONUNIT_DIAG, A, B, IERROR )
C
C Stop timing 
C
         CALL MPI_BARRIER( MPI_COMM_WORLD, IERROR )
         TIME = MPI_WTIME() - TIME
C
C Report performance
C
         IF ( ME .EQ. 0 ) THEN
            CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, IERROR )

            OPERATION_COUNT = 2.0/3.0 * SIZE * SIZE * SIZE

            IF ( MPI_DOUBLE_COMPLEX .EQ. DATATYPE ) 
     &           OPERATION_COUNT = 4 * OPERATION_COUNT

            PRINT *, "N = ", SIZE, " TIME = ", TIME, " MFLOPS/node = ", 
     &           OPERATION_COUNT / TIME * 1.0e-6 / NPROCS
         ENDIF
C
C    /* Process the answer.  As an example, this routine brings 
C       result x (stored in b) to processor 0 and prints first and 
C       last entry */
C
         CALL PROCESS_ANSWER_F( B )
C
C    /* Check answer by overwriting b <- b - x (where b holds computed
C       approximation to x) */
C     
         CALL PLA_Axpy_F( MINUS_ONE, X, B, IERROR )
C
C    /* Create 1x1 multiscalars to hold largest (in abs. value) element 
C       of b - x and index of largest value */
C
         CALL PLA_MscalarCreate_F( DATATYPE, 
     &                              PLA_ALL_ROWS, PLA_ALL_COLS,
     &                              1, 1, TEMPL, B_NORM, IERROR )
C
         CALL PLA_MscalarCreate_F( MPI_INTEGER, 
     &                              PLA_ALL_ROWS, PLA_ALL_COLS,
     &                              1, 1, TEMPL, INDEX, IERROR )
C
         CALL PLA_Iamax_F ( B, INDEX, B_NORM, IERROR )
C
C    /* Report norm of b - x */
C
         IF ( ME .EQ. 0 ) THEN
            CALL PLA_ObjGetLocalContents_F( 
     &           B_NORM, PLA_NO_TRANS, DUMMY, DUMMY,
     &           B_NORM_VALUE, 1, 1, IERROR )
            PRINT *, "Largest element in (computed x) - x : ", 
     &               B_NORM_VALUE
         ENDIF
      ENDIF
C
C  /* Free the linear algebra objects */
C
      CALL PLA_ObjFree_F( A, IERROR )
      CALL PLA_ObjFree_F( X, IERROR )
      CALL PLA_ObjFree_F( B, IERROR )
      CALL PLA_ObjFree_F( B_NORM, IERROR )
      CALL PLA_ObjFree_F( PIVOTS, IERROR )
      CALL PLA_ObjFree_F( INDEX, IERROR )
      CALL PLA_ObjFree_F( MINUS_ONE, IERROR )
      CALL PLA_ObjFree_F( ZERO, IERROR )
      CALL PLA_ObjFree_F( ONE, IERROR )
C
C  /* Free the template */
C
      CALL PLA_TempFree_F( TEMPL, IERROR )
C
C  /* Finalize PLAPACK and MPI */
C
      CALL PLA_Finalize_F(IERROR )
      CALL MPI_Finalize( IERROR )

      STOP
      END




