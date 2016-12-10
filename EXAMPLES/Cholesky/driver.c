/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

#if MANUFACTURE == CRAY
#include <mpp/rastream.h>
#endif  

#define datatype MPI_DOUBLE_COMPLEX

int PLA_Chol( int, PLA_Obj);

int main(int argc, char *argv[])
{
  MPI_Comm comm;
  PLA_Template template = NULL;
  PLA_Obj  A = NULL,           x = NULL,         b = NULL, 
           minus_one = NULL,   b_norm = NULL;
  double   *locA, b_norm_value, time, flops;
  int      size, nb_distr, nb_alg, me, nprocs, dummy, info = 0;
  int      nprows, npcols, ierror;

  PLA_Obj A11 = NULL;

  MPI_Init(&argc, &argv);

#if MANUFACTURE == CRAY
  set_d_stream( 1 );
#endif

  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


  if (0 == me) {
    printf("enter nprows, npcols:\n");
    scanf("%d %d", &nprows, &npcols );
    printf("enter matrix size, distr. block size, and alg. block size:\n");
    scanf("%d %d %d", &size, &nb_distr, &nb_alg );
    printf("Turn on error checking? (1 = YES, 0 = NO):\n");
    scanf("%d", &ierror );
  }

  MPI_Bcast(&nprows,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&npcols,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size,     1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_alg,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ierror, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if ( ierror ) 
    PLA_Set_error_checking( ierror, TRUE, TRUE, FALSE );
  else
    PLA_Set_error_checking( ierror, FALSE, FALSE, FALSE );

  pla_Environ_set_nb_alg (PLA_OP_ALL_ALG,
			  nb_alg);

/*  PLA_Comm_1D_to_2D_ratio(MPI_COMM_WORLD, 1.0, &comm); */
  PLA_Comm_1D_to_2D(MPI_COMM_WORLD, nprows, npcols, &comm);
  PLA_Init(comm);

  PLA_Temp_create( nb_distr, 0, &template );

  PLA_Matrix_create(  datatype, size, size, 
                      template, PLA_ALIGN_FIRST, PLA_ALIGN_FIRST, &A );

  PLA_Mvector_create( datatype, size, 1, template, PLA_ALIGN_FIRST, &x );

  PLA_Mvector_create( datatype, size, 1, template, PLA_ALIGN_FIRST, &b );

  /* Create usual duplicated scalar constants with same datatype and
     template as A_panel */
  PLA_Create_constants_conf_to( A, &minus_one, NULL, NULL );

  PLA_Mscalar_create( MPI_DOUBLE, PLA_ALL_ROWS, PLA_ALL_COLS,
		      1, 1, template, &b_norm );

  create_problem( A, x, b ); 

  pla_Environ_set_nb_alg( PLA_OP_ALL_ALG, nb_alg);

  MPI_Barrier( MPI_COMM_WORLD );
  time = MPI_Wtime( );

  info = PLA_Chol ( PLA_LOWER_TRIANGULAR, A );

  if ( info != 0 ) {
    printf("nonpositive definite matrix detected at row %d\n", info);
  }
  else {
    PLA_Trsv( PLA_LOWER_TRIANGULAR, PLA_NO_TRANSPOSE, PLA_NONUNIT_DIAG, A, b ); 

    PLA_Trsv( PLA_LOWER_TRIANGULAR, 
	       ( datatype == MPI_DOUBLE || datatype == MPI_FLOAT ?
		 PLA_TRANSPOSE : PLA_CONJUGATE_TRANSPOSE ),    
	       PLA_NONUNIT_DIAG, A, b );

    time = MPI_Wtime() - time;

    if ( me == 0 )  {
      flops = 1.0/3.0 * size * size * size;
      if ( datatype == MPI_COMPLEX || datatype == MPI_DOUBLE_COMPLEX )
	flops *= 4;

      printf("n = %d, time = %lf, MFLOPS/node = %lf\n", size, time,
	     flops / time * 1.0e-6 / nprocs );
    }

    PLA_Axpy( minus_one, x, b );

    PLA_Nrm2( b, b_norm );

    if ( me == 0 ) {
      PLA_Obj_get_local_contents( b_norm, PLA_NO_TRANS, &dummy, &dummy,
				  &b_norm_value, 1, 1 );
      printf( "Norm of difference between original x and computed x: %le\n", 
	      b_norm_value );
    }
  } 

  PLA_Obj_free(&A);            PLA_Obj_free(&x);
  PLA_Obj_free(&b);            PLA_Obj_free(&minus_one);
  PLA_Obj_free(&b_norm);

  PLA_Temp_free(&template);
  PLA_Finalize( );
  MPI_Finalize( );
}
