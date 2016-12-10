      SUBROUTINE CREATE_PROBLEM_F( A, X, B )
C
C /* Compute random A and x and let b = A x. */
C
C
      IMPLICIT NONE

      INCLUDE "mpif.h"
      INCLUDE "PLAf.h"
C
C PLA_Obj Parameters
C
      INTEGER 
     &     A,
     &     X,
     &     B
C
C PLA_Obj Local
C
      INTEGER
     &     MINUS_ONE,
     &     ZERO,
     &     ONE
C
C Maximum matrix dimension (length of column)
C
      INTEGER MAX_SIZE
      PARAMETER ( MAX_SIZE = 10000 )
C
      INTEGER SIZE, ME, NPROCS, I, J, IERROR
      DOUBLE PRECISION Z_COL( MAX_SIZE )
      DOUBLE PRECISION D_ONE
C
C MPI_Datatype
C
      INTEGER DATATYPE
C
C Initialize objects to NULL
C
      MINUS_ONE = 0
      ZERO = 0
      ONE  = 0
      
      D_ONE = 1.0D00
C
C /* Extract dimension of A */
C
      CALL PLA_ObjGlobalLength_F( A, SIZE, IERROR )
      IF ( SIZE .GT. MAX_SIZE ) THEN
         PRINT *, "SIZE > MAX_SIZE IN CREATE_PROBLEM_F"
         STOP
      ENDIF

      CALL MPI_COMM_RANK( MPI_COMM_WORLD, ME, IERROR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NPROCS, IERROR )

      CALL PLA_ObjSetToZero_F( A, IERROR )
      CALL PLA_ObjSetToZero_F( X, IERROR )
      CALL PLA_ObjSetToZero_F( B, IERROR )

      CALL PLA_ObjDatatype_F( A, DATATYPE, IERROR )
C
C /* Enter Application Interface mode */
C
      CALL PLA_ApiBegin_F( IERROR )
C
C /* Open objects A and x for read/write */
C
         CALL PLA_ObjApiOpen_F( A, IERROR )
         CALL PLA_ObjApiOpen_F( X, IERROR )
C
C /* Seed the random number generator */
C
         CALL SEED_RANDOM( ME * 1793 )
C
C /* Columns are computed by processors in a round-robin fashion */
C
         DO J=ME+1, SIZE, NPROCS
C
C      /* Fill column j with random numbers */
C
            DO I=1, SIZE
               IF ( MPI_DOUBLE_COMPLEX .EQ. DATATYPE ) THEN
                  CALL COMPLEX_RANDOM( Z_COL( I ) )
               ELSE IF ( MPI_DOUBLE_PRECISION .EQ. DATATYPE ) THEN
                  CALL REAL_RANDOM( Z_COL( I ) )
               ELSE
                  PRINT *, "DATATYPE UNKNOWN"
                  STOP
               ENDIF
            ENDDO
C
C      /* Submit the column to the global matrix */
C
            CALL PLA_ApiAxpyMatrixToGlobal_F(
     &           SIZE, 1, D_ONE, Z_COL, SIZE, A, 0, J-1, IERROR )
         ENDDO
         IF ( ME .EQ. 0 ) THEN
            DO I=1, SIZE
               IF ( MPI_DOUBLE_COMPLEX .EQ. DATATYPE ) THEN
                  CALL COMPLEX_RANDOM( Z_COL( I ) )
               ELSE IF ( MPI_DOUBLE_PRECISION .EQ. DATATYPE ) THEN
                  CALL REAL_RANDOM( Z_COL( I ) )
               ELSE
                  PRINT *, "DATATYPE UNKNOWN"
                  STOP
               ENDIF	
            ENDDO
            CALL PLA_ApiAxpyVectorToGlobal_F( 
     &           SIZE, D_ONE, Z_COL( 1 ), 1, X, 0, IERROR )
         ENDIF
C         
C    /* Close the objects */
C
      CALL PLA_ObjApiClose_F( X, IERROR )
      CALL PLA_ObjApiClose_F( A, IERROR )
C
C                         /* leave Application Interface mode */
C
      CALL PLA_ApiEnd_F( IERROR )
C
C  /* Compute b = A x by performing a matrix-vector multiplication */
C
      CALL PLA_CreateConstantsConfTo_F( A, MINUS_ONE, ZERO, ONE, 
     &     IERROR )
      CALL PLA_Gemv_F( PLA_NO_TRANSPOSE, ONE, A, X, ZERO, B, IERROR )
C
C  /* Free the temporary objects */
C
      CALL PLA_ObjFree_F( MINUS_ONE, IERROR )
      CALL PLA_ObjFree_F( ZERO, IERROR )
      CALL PLA_ObjFree_F( ONE, IERROR )

      RETURN
      END
      
