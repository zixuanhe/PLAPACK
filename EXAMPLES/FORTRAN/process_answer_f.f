      SUBROUTINE PROCESS_ANSWER_F( X )
      
      IMPLICIT NONE
      
      INCLUDE "mpif.h"
      INCLUDE "PLAf.h"
C 
C PLA_OBJ
C
      INTEGER X
C
      INTEGER SIZE, ME, I, IERROR
C
      DOUBLE PRECISION D_ONE
      PARAMETER ( D_ONE = 1.0 )
C
      INTEGER MAX_SIZE
      PARAMETER ( MAX_SIZE = 10000 )
C
      DOUBLE PRECISION SOLUTION( MAX_SIZE )
C
C MPI Datatype
C
      INTEGER DATATYPE

      CALL PLA_ObjDatatype_F( X, DATATYPE, IERROR )
      CALL PLA_ObjGlobalLength_F( X, SIZE, IERROR )
      IF ( SIZE .GT. MAX_SIZE ) THEN
         PRINT *, "SIZE too large in process_answer_f"
         STOP
      ENDIF
C
C  /* What is my rank in the communicator? */
C
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, ME, IERROR )
C
C  /* Begin Application Interface Mode */
C
      CALL PLA_ApiBegin_F( IERROR )
C
C    /* Open object x for read/write mode */
C
         CALL PLA_ObjApiOpen_F( X, IERROR )
         
            IF ( ME .EQ. 0 ) THEN
C
C Initialize solution to zero
C           
               IF ( DATATYPE .EQ. MPI_COMPLEX ) THEN
                  DO I=1, SIZE
                     SOLUTION( I ) = ( 0.0E00, 0.0E00 )
                  ENDDO
               ELSE IF ( DATATYPE .EQ. MPI_DOUBLE_PRECISION ) THEN
                  DO I=1, SIZE
                     SOLUTION( I ) = 0.0E00
                  ENDDO
               ELSE
                  PRINT *, "Datatype unknown in process_answer"
               ENDIF
C
C      /* Initiate extraction of contents of x to local buffer */
C
               CALL PLA_ApiAxpyGlobalToVector_F( 
     &              SIZE, D_ONE, X, 0, SOLUTION, 1, IERROR )
            ENDIF
            
C  
C  /* Wait until data actually arrives */
C
            CALL PLA_ObjApiSync_F( X, IERROR )
            
            IF ( ME .EQ. 0 ) THEN
               PRINT *, "First computed entry:", SOLUTION( 1 )
               PRINT *, "Last  computed entry:", SOLUTION( SIZE )
            ENDIF
C
C    /* Close object x */
C
         CALL PLA_ObjApiClose_F( X, IERROR )
C
C  /* Exit API mode */
C
      CALL PLA_ApiEnd_F( IERROR )

      RETURN
      END
