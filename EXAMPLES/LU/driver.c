/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/*
 *   PLAPACK General Linear Solver example application.
 *
 *   Purpose:   Implement the solution of a problem requiring
 *              a general linear solver, including the problem
 *              generation phase, the solution of the system,
 *              and the retrieval of the solution data.
 *
 *   Datatypes: MPI_FLOAT, MPI_DOUBLE, MPI_DOUBLE_COMPLEX 
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
    A = NULL,           
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

  /* Two different ways to create a 2D communicator */
/*  PLA_Comm_1D_to_2D_ratio(MPI_COMM_WORLD, 0.25, &comm); */ 
  PLA_Comm_1D_to_2D(MPI_COMM_WORLD, nprows, npcols, &comm); 

  /* Initialize PLAPACK */
  PLA_Init(comm);

  /* Create object distribution template */
  PLA_Temp_create( nb_distr, 0, &templ );

  /* Set the datatype : MPI_DOUBLE, MPI_FLOAT or MPI_DOUBLE_COMPLEX */
  datatype = MPI_DOUBLE_COMPLEX;

  /* Create objects for problem to be solved */
  PLA_Matrix_create(  datatype, size, size, 
                      templ, PLA_ALIGN_FIRST, PLA_ALIGN_FIRST, &A );

  PLA_Mvector_create( datatype, size, 1, templ, PLA_ALIGN_FIRST, &x );

  PLA_Mvector_create( datatype, size, 1, templ, PLA_ALIGN_FIRST, &b );
   
  PLA_Mvector_create( MPI_INT, size, 1, templ, PLA_ALIGN_FIRST, &pivots );
   
  /* Create 1x1 multiscalars to hold largest (in abs. value) element 
     of b - x and index of largest value */
  if ( datatype == MPI_DOUBLE_COMPLEX || datatype == MPI_DOUBLE )
    PLA_Mscalar_create( MPI_DOUBLE,
			 PLA_ALL_ROWS, PLA_ALL_COLS,
			 1, 1, templ, &b_norm );
  else   if ( datatype == MPI_COMPLEX || datatype == MPI_FLOAT )
    PLA_Mscalar_create( MPI_FLOAT,
			 PLA_ALL_ROWS, PLA_ALL_COLS,
			 1, 1, templ, &b_norm );

  /* Create duplicated scalar constants with same datatype and template as A */
  PLA_Create_constants_conf_to( A, &minus_one, NULL, NULL );

  /* Create a problem to be solved: A x = b */
  create_problem( A, x, b );

  /* Start timing */
  MPI_Barrier( MPI_COMM_WORLD );
  time = MPI_Wtime( );

  /* Factor P A -> L U overwriting lower triangular portion of A with L, upper, U */
  info = PLA_LU ( A, pivots);

  if ( info != 0 ) {
    printf("Zero pivot encountered at row %d.\n", info);
  }
  else {
    /* Apply the permutations to the right hand sides */
    PLA_Apply_pivots_to_rows ( b, pivots);

    /* Solve L y = b, overwriting b with y */
    PLA_Trsv( PLA_LOWER_TRIANGULAR, PLA_NO_TRANSPOSE, PLA_UNIT_DIAG, A, b );

    /* Solve U x = y (=b), overwriting b with x */
    PLA_Trsv( PLA_UPPER_TRIANGULAR, PLA_NO_TRANSPOSE,  PLA_NONUNIT_DIAG, A, b );

    /* Stop timing */
    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime() - time;

    /* Report performance */
    if ( me == 0 ) {
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

      operation_count = 2.0/3.0 * size * size * size;

      if ( MPI_DOUBLE_COMPLEX == datatype )
	operation_count *= 4;

      switch( datatype ){
      case MPI_DOUBLE_COMPLEX: 
	printf( "datatype timed: double precision complex\n"); 
	printf( "(op count = 8/3 n^3)\n"); 
	break;
      case MPI_DOUBLE: 
	printf( "datatype timed: double precision\n"); 
	printf( "(op count = 2/3 n^3)\n"); 
	break;
      case MPI_COMPLEX: 
	printf( "datatype timed: single precision complex\n"); 
	printf( "(op count = 8/3 n^3)\n"); 
	break;
      case MPI_FLOAT: 
	printf( "datatype timed: single precision\n"); 
	printf( "(op count = 2/3 n^3)\n"); 
	break;
      }

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
  PLA_Obj_free(&b_norm);       PLA_Obj_free (&pivots);

  /* Free the template */
  PLA_Temp_free(&templ);

  /* Finalize PLAPACK and MPI */
  PLA_Finalize( );
  MPI_Finalize( );
}




