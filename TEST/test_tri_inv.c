/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "mpi.h"
#include "PLA.h"

int main(int argc, char *argv[])
{
  MPI_Comm 
    comm = MPI_COMM_NULL;
  MPI_Datatype
    datatype;
  PLA_Template 
    templ = NULL;
  PLA_Obj  
    A      = NULL, A_orig = NULL, residual = NULL,
    minus_one   = NULL, one  = NULL, diff = NULL;
  int      
    n,
    nb_distr, nb_alg,
    error, parameters, sequential,
    me, nprocs, nprows, npcols,
    itype;
  double 
    time,
    flops;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (me==0) {
    printf("enter mesh size:\n");
    scanf("%d%d", &nprows, &npcols );
    printf("mesh size  = %d x %d \n", nprows, npcols );
    printf("enter distr. block size:\n");
    scanf("%d", &nb_distr );
    printf("nb_distr = %d\n", nb_distr );
    printf("enter alg. block size:\n");
    scanf("%d", &nb_alg );
    printf("nb_alg = %d\n", nb_alg );
    printf("turn on error checking? (0 = NO, 1 = YES):\n");
    scanf("%d", &error );
    printf("error checking = %d\n", error );
    printf("turn on parameter checking? (0 = NO, 1 = YES):\n");
    scanf("%d", &parameters );
    printf("parameter checking = %d\n", parameters );
    printf("turn on sequential checking? (0 = NO, 1 = YES):\n");
    scanf("%d", &sequential );
    printf("sequential checking = %d\n", sequential );
  }
  MPI_Bcast(&nprows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&npcols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_alg,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&parameters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sequential, 1, MPI_INT, 0, MPI_COMM_WORLD);

  pla_Environ_set_nb_alg( PLA_OP_ALL_ALG, nb_alg );

  PLA_Set_error_checking( error, parameters, sequential, FALSE );

/*  PLA_Comm_1D_to_2D_ratio(MPI_COMM_WORLD, 1.0, &comm); */
  PLA_Comm_1D_to_2D(MPI_COMM_WORLD, nprows, npcols, &comm); 

  PLA_Init(comm);
    
  PLA_Temp_create( nb_distr, 0, &templ );
    
  while ( TRUE ){
    if (me==0) {
      printf("enter datatype:\n");
      printf("-1 = quit\n");
      printf(" 0 = float\n");
      printf(" 1 = double\n");
      printf(" 2 = complex\n");
      printf(" 3 = double complex\n");
      scanf("%d", &itype );
      printf("itype = %d\n", itype );
    }
    MPI_Bcast(&itype, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if ( itype == -1 ) break;
    switch( itype ){
    case 0:
      datatype = MPI_FLOAT;
      break;
    case 1:
      datatype = MPI_DOUBLE;
      break;
    case 2:
      datatype = MPI_COMPLEX;
      break;
    case 3:
      datatype = MPI_DOUBLE_COMPLEX;
      break;
    default:
      PLA_Abort( "unknown datatype", __LINE__, __FILE__ );
    }

    if (me==0) {
      printf("enter n:\n");
      scanf("%d", &n );
      printf("n = %d\n", n );
    }

    MPI_Bcast(&n,     1, MPI_INT, 0, MPI_COMM_WORLD);

    PLA_Matrix_create( datatype, 
		        n, 
		        n,
			templ, 
                        PLA_ALIGN_FIRST, 
                        PLA_ALIGN_FIRST, 
                        &A ); 

    PLA_Matrix_create_conf_to( A, &A_orig );

    PLA_Matrix_create_conf_to( A, &residual );

    create_problem( A );

    {
      double d_n;

      d_n = (double) n;
      PLA_Shift( A, MPI_DOUBLE, &d_n );
    }

    PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

    PLA_Local_copy( A, A_orig );

    /* Use invert routine that uses factors */
    
    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime ();

    PLA_Triangular_invert( PLA_LOWER_TRIANGULAR, A );

    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime () - time;

    flops = 2.0/3.0 * n * n * n;

    if ( me == 0 ) 
      printf("%d time = %f, MFLOPS/node = %10.4lf \n", n, time,
	       flops / time * 1.0e-6 / nprocs );

    PLA_Obj_set_to_identity( residual );

    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, 
			     A_orig );
    PLA_Set_triang_to_zero( PLA_LOWER_TRIANGULAR, PLA_NONUNIT_DIAG, 
			     A );

    PLA_Gemm( PLA_NO_TRANS, PLA_NO_TRANS, minus_one, A_orig, A, one, 
	       residual );

    PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
			 1, 1, templ, &diff );

    PLA_Matrix_one_norm( residual, diff );
  }

  PLA_Obj_free( &A );
  PLA_Obj_free( &A_orig );
  PLA_Obj_free( &residual );
  PLA_Obj_free( &minus_one );
  PLA_Obj_free( &one );
  PLA_Obj_free( &diff );

  PLA_Temp_free(&templ);
  PLA_Finalize( );
  MPI_Finalize( );
  
}
