/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/*
 *   PLAPACK LINPACK benchmark driver
 *
 */

#include "PLA.h"                 /* PLAPACK header file */   

#if MANUFACTURE == CRAY
#include <mpp/rastream.h>
#endif  

int main(int argc, char *argv[])
{
  /* Declarations */
  MPI_Comm 
    comm;

  PLA_Template 
    templ = NULL;

  PLA_Obj  
    A = NULL,    rhs = NULL,    
    A_append = NULL,
    pivots = NULL,
    x = NULL,         
    b = NULL, 
    b_norm = NULL,
    index = NULL,
    minus_one = NULL;

  double 
    operation_count,
    b_norm_value, 
    time;

  int  
    size, 
    nb_distr, nb_alg, 
    me, nprocs, 
    nprows, npcols,
    dummy, 
    ierror,
    info = 0;
  
  MPI_Datatype 
    datatype;

  /* Initialize MPI */
  MPI_Init(&argc, &argv);

#if MANUFACTURE == CRAY
  set_d_stream( 1 );
#endif

  /* Get problem size and distribution block size and broadcast */
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  if (0 == me) {
    printf("enter processor mesh dimension ( rows cols ):\n");
    scanf("%d %d", &nprows, &npcols );
    printf("enter matrix size, distr. block size:\n");
    scanf("%d %d", &size, &nb_distr );
    printf("enter algorithmic blocksize:\n");
    scanf("%d", &nb_alg );
    printf("Turn on error checking? (1 = YES, 0 = NO):\n");
    scanf("%d", &ierror );
  }

  MPI_Bcast(&nprows,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&npcols,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size,     1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_alg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ierror, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if ( ierror ) 
    PLA_Set_error_checking( ierror, TRUE, TRUE, FALSE );
  else
    PLA_Set_error_checking( ierror, FALSE, FALSE, FALSE );

  pla_Environ_set_nb_alg (PLA_OP_ALL_ALG,
			  nb_alg);

  /* Create a 2D communicator */
  PLA_Comm_1D_to_2D(MPI_COMM_WORLD, nprows, npcols, &comm); 

  /* Initialize PLAPACK */
  PLA_Init(comm);

  /* Create object distribution template */
  PLA_Temp_create( nb_distr, 0, &templ );

  /* Set the datatype */
  datatype = MPI_DOUBLE;

  /* Create objects for problem to be solved */

  /* Matrix A is big enough to hold the right-hand-side appended */
  PLA_Matrix_create(  datatype, size, size+1, 
		      templ, PLA_ALIGN_FIRST, PLA_ALIGN_FIRST, &A_append );

  PLA_Mvector_create( datatype, size, 1, templ, PLA_ALIGN_FIRST, &x );

  PLA_Mvector_create( datatype, size, 1, templ, PLA_ALIGN_FIRST, &b );
   
  PLA_Mvector_create( MPI_INT, size, 1, templ, PLA_ALIGN_FIRST, &pivots );
   
  /* Create 1x1 multiscalars to hold largest (in abs. value) element 
     of b - x and index of largest value */
  PLA_Mscalar_create( MPI_DOUBLE,
		      PLA_ALL_ROWS, PLA_ALL_COLS,
			 1, 1, templ, &b_norm );

  /* Create duplicated scalar constants with same datatype and template as A */
  PLA_Create_constants_conf_to( A_append, &minus_one, NULL, NULL );

  /* View the appended system as the matrix and the right-hand-side */
  PLA_Obj_vert_split_2( A_append, -1, &A, &rhs );

  /* Create a problem to be solved: A x = b */
  create_problem( A, x, b );

  /* Copy b to the appended column */
  PLA_Copy( b, rhs );

  /* Start timing */
  MPI_Barrier( MPI_COMM_WORLD );
  time = MPI_Wtime( );

  /* Factor P A_append -> L U overwriting lower triangular portion of A with L, upper, U */

  info = PLA_LU( A_append, pivots);

  if ( info != 0 ) {
    printf("Zero pivot encountered at row %d.\n", info);
  }
  else {
    /* Apply the permutations to the right hand sides */
    /* Not necessery since system was appended */
    /* PLA_Apply_pivots_to_rows ( b, pivots); */
    
    /* Solve L y = b, overwriting b with y */
    /* Not necessary since the system was appended */
    /* PLA_Trsv( PLA_LOWER_TRIANGULAR, PLA_NO_TRANSPOSE, PLA_UNIT_DIAG, A, b ); */
    PLA_Copy( rhs, b );

    /* Solve U x = y (=b), overwriting b with x */
    PLA_Trsv( PLA_UPPER_TRIANGULAR, PLA_NO_TRANSPOSE,  PLA_NONUNIT_DIAG, A, b );

    /* Stop timing */
    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime() - time;

    /* Report performance */
    if ( me == 0 ) {
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

      operation_count = 2.0/3.0 * size * size * size;

      printf("n = %d, time = %lf, MFLOPS/node = %lf\n", size, time,
	      operation_count / time * 1.0e-6 / nprocs );
    }

    /* Process the answer.  As an example, this routine brings 
       result x (stored in b) to processor 0 and prints first and 
       last entry */
    Process_answer( b ); 

    /* Check answer by overwriting b <- b - x (where b holds computed
       approximation to x) */

    PLA_Axpy( minus_one, x, b );

    PLA_Nrm2( b, b_norm);

    /* Report norm of b - x */
    if ( me == 0 ) {
      PLA_Obj_get_local_contents( b_norm, PLA_NO_TRANS, &dummy, &dummy,
				  &b_norm_value, 1, 1 );
      printf( "Norm2 of x - computed x : %le\n", b_norm_value );
    }
  } 

  /* Free the linear algebra objects */
  PLA_Obj_free(&A);            PLA_Obj_free(&x);
  PLA_Obj_free(&b);            PLA_Obj_free(&minus_one);
  PLA_Obj_free(&b_norm);       PLA_Obj_free(&pivots);
  PLA_Obj_free(&A_append);     PLA_Obj_free(&rhs);

  /* Free the template */
  PLA_Temp_free(&templ);

  /* Finalize PLAPACK and MPI */
  PLA_Finalize( );
  MPI_Finalize( );
}




