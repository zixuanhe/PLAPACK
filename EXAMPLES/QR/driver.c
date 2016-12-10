/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/*
 *   PLAPACK Least squares solver example application
 *
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
    Q = NULL,
    s = NULL,
    pivots = NULL,
    x = NULL,         
    b = NULL, 
    b_norm = NULL,
    one = NULL,
    minus_one = NULL;

  double 
    operation_count,
    b_norm_value, 
    time;

  int  
    size, 
    nb_distr, 
    nb_alg,
    me, nprocs, 
    nprows, npcols,
    dummy, 
    info = 0,
    form_Q,
    ierror;
  
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
    printf("enter matrix size, distribution block size, algorithmic block size:\n");
    scanf("%d %d %d", &size, &nb_distr, &nb_alg );
    printf("If Q is to be formed, enter 1.  O/w enter 0:\n");
    scanf("%d", &form_Q );
    printf("Turn on error checking? (1 = YES, 0 = NO):\n");
    scanf("%d", &ierror );
  }

  MPI_Bcast(&nprows,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&npcols,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size,     1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_alg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&form_Q, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ierror, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if ( ierror ) 
    PLA_Set_error_checking( ierror, TRUE, TRUE, FALSE );
  else
    PLA_Set_error_checking( ierror, FALSE, FALSE, FALSE );

  /* Two different ways to create a 2D communicator */
/*  PLA_Comm_1D_to_2D_ratio(MPI_COMM_WORLD, 0.25, &comm); */ 
  PLA_Comm_1D_to_2D(MPI_COMM_WORLD, nprows, npcols, &comm); 

  /* Initialize PLAPACK */
  PLA_Init(comm);

/*  PLA_Set_check_parameters ( PLA_CHECK_PARAMETER ); */

  /* Create object distribution template */
  PLA_Temp_create( nb_distr, 0, &templ );

  pla_Environ_set_nb_alg (PLA_OP_ALL_ALG,
			  nb_alg);

  /* Set the datatype : MPI_DOUBLE, MPI_FLOAT or MPI_DOUBLE_COMPLEX */
  datatype = MPI_DOUBLE;

  /* Create objects for problem to be solved */
  PLA_Matrix_create(  datatype, size, size, 
                      templ, PLA_ALIGN_FIRST, PLA_ALIGN_FIRST, &A ); 
  PLA_Mvector_create(  datatype, size, 1,
                      templ, PLA_ALIGN_FIRST, &s );

  if ( form_Q )
    PLA_Matrix_create(  datatype, size, size, 
                        templ, PLA_ALIGN_FIRST, PLA_ALIGN_FIRST, &Q );

  PLA_Mvector_create( datatype, size, 1, templ, PLA_ALIGN_FIRST, &x );

  PLA_Mvector_create( datatype, size, 1, templ, PLA_ALIGN_FIRST, &b );
   
  /* Create duplicated scalar constants with same datatype and template as A */
  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  /* Create a problem to be solved: A x = b */
  create_problem( A, x, b );

  /* Start timing */
  MPI_Barrier( MPI_COMM_WORLD );
  time = MPI_Wtime( );

  /* Factor P A -> L U overwriting lower triangular portion of A with L, upper, U */
  info = PLA_Least_sqr_solve ( A, s, ( form_Q ? Q : PLA_DUMMY ), b); 

  time = MPI_Wtime() - time;

/*  PLA_Trmm( PLA_SIDE_RIGHT, PLA_UPPER_TRIANGULAR, PLA_NO_TRANSPOSE, 
            PLA_NONUNIT_DIAG, one, A, Q );  */

  /* Report performance */
  if ( me == 0 ) {
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    operation_count = 4.0/3.0 * size * size * size;
    if ( form_Q ) operation_count *= 2;

    printf("n = %d, time = %lf, MFLOPS/node = %lf\n", size, time,
	      operation_count / time * 1.0e-6 / nprocs );
    }

  /* Check answer by overwriting b <- b - x (where b holds computed
     approximation to x) */
  PLA_Axpy( minus_one, x, b ); 

  PLA_Mscalar_create((datatype==MPI_DOUBLE_COMPLEX?MPI_DOUBLE:datatype), 
		     PLA_ALL_ROWS, PLA_ALL_COLS,
		     1, 1, templ, &b_norm );

  PLA_Nrm2( b, b_norm);

  /* Report norm of b - x */
  if ( me == 0 ) {
    PLA_Obj_get_local_contents( b_norm, PLA_NO_TRANS, &dummy, &dummy,
			       &b_norm_value, 1, 1 );
    printf( "Norm2 of x - computed x : %le\n", b_norm_value );
    }

  /* Free the linear algebra objects */
  PLA_Obj_free(&A);            PLA_Obj_free(&x);
  PLA_Obj_free(&b);            PLA_Obj_free(&minus_one);
  PLA_Obj_free(&b_norm);       PLA_Obj_free (&one);
  PLA_Obj_free (&Q);
  PLA_Obj_free( &s );

  /* Free the template */
  PLA_Temp_free(&templ);

  /* Finalize PLAPACK and MPI */
  PLA_Finalize( );
  MPI_Finalize( );
}




