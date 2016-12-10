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
    A      = NULL, B    = NULL, C = NULL, 
    alpha  = NULL, beta = NULL;
  int      
    m, n, k, 
    nb_distr, nb_alg,
    error, parameters, sequential,
    me, nprocs, 
    itype,
    itransa, itransb, transa, transb,
    variant, version; 
  double 
    time,
    flops;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (me==0) {
    printf("enter distr. block size:\n");
    scanf("%d", &nb_distr );
    printf("nb_distr = %d\n", nb_distr );
    printf("enter algorithmic block size:\n");
    scanf("%d", &nb_alg );
    printf("nb_alg = %d\n", nb_distr );
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
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_alg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&parameters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sequential, 1, MPI_INT, 0, MPI_COMM_WORLD);

  pla_Environ_set_nb_alg( PLA_OP_ALL_ALG, nb_alg );

  PLA_Set_error_checking( error, parameters, sequential, FALSE);

  PLA_Comm_1D_to_2D_ratio(MPI_COMM_WORLD, 1.0, &comm);
  PLA_Init( comm );
    
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
      printf("variant? (0 = A, 1 = B, 2 = C, 3 = standard):\n");
      scanf("%d", &variant );
      printf("variant = %d\n", variant );
      printf("enter m, n, k:\n");
      scanf("%d %d %d", &m, &n, &k );
      printf("m = %d, n = %d, k = %d\n", m , n, k );
      printf("enter transa, transb (0=notrans, 1=trans, 2=conj trans):\n");
      scanf("%d %d", &itransa, &itransb );
      printf("itransa = %d, itransb = %d\n", itransa , itransb );
    }

    MPI_Bcast(&variant, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&itransa,  1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&itransb,  1, MPI_INT, 0, MPI_COMM_WORLD);

    switch( itransa ){
    case 0: 
      transa = PLA_NO_TRANS;
      break;
    case 1: 
      transa = PLA_TRANS;
      break;
    case 2: 
      transa = PLA_CONJ_TRANS;
      break;
    default:
      PLA_Abort( "transpose case not supported", __LINE__, __FILE__ );
    }

    switch( itransb ){
    case 0: 
      transb = PLA_NO_TRANS;
      break;
    case 1: 
      transb = PLA_TRANS;
      break;
    case 2: 
      transb = PLA_CONJ_TRANS;
      break;
    default:
      PLA_Abort( "transpose case not supported", __LINE__, __FILE__ );
    }

    switch( variant ){
    case 0: 
      version = PLA_VERSION_A;
      break;
    case 1: 
      version = PLA_VERSION_B;
      break;
    case 2: 
      version = PLA_VERSION_C;
      break; 
    case 3:
      break;
    default:
      PLA_Abort( "variant not supported", __LINE__, __FILE__ );
    }

    PLA_Matrix_create( datatype, 
		       ( transa == PLA_NO_TRANS ? m : k ),
		       ( transa == PLA_NO_TRANS ? k : m ),
			templ, 
                        PLA_ALIGN_FIRST, 
                        PLA_ALIGN_FIRST, 
                        &A ); 

    create_problem( A );

    PLA_Matrix_create( datatype, 
		       ( transb == PLA_NO_TRANS ? k : n ),
		       ( transb == PLA_NO_TRANS ? n : k ),
		        templ, 
                        PLA_ALIGN_FIRST, 
                        PLA_ALIGN_FIRST, 
                        &B ); 
      
    create_problem( B );

    PLA_Matrix_create( datatype, 
		        m, 
                        n,
			templ, 
                        PLA_ALIGN_FIRST, 
                        PLA_ALIGN_FIRST, 
                        &C ); 

    create_problem( C );

    PLA_Mscalar_create( datatype, 
		        PLA_ALL_ROWS,
                        PLA_ALL_COLS,
			1, 
			1,
			templ, 
                        &alpha ); 

    PLA_Obj_set_to_one( alpha );

    PLA_Mscalar_create( datatype, 
		        PLA_ALL_ROWS,
                        PLA_ALL_COLS,
			1, 
			1,
			templ, 
                        &beta ); 

    PLA_Obj_set_to_zero( beta );

    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime ();

    if ( variant == 3 )
      PLA_Gemm( transa, transb, alpha, A, B, beta, C );
    else
      PLA_Gemm_x( version, transa, transb, alpha, A, B, beta, C );

    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime () - time;

    flops = 2.0 * m * n * k;
    if ( datatype == MPI_COMPLEX || datatype == MPI_DOUBLE_COMPLEX )
      flops *= 4;

    if ( me == 0 ) 
      printf("%d %d %d time = %f, MFLOPS/node = %10.4lf \n", m, n, k, time,
	       flops / time * 1.0e-6 / nprocs );
  }

  PLA_Obj_free( &A );
  PLA_Obj_free( &B );
  PLA_Obj_free( &C );
  PLA_Obj_free(&alpha);
  PLA_Obj_free(&beta);

  PLA_Temp_free(&templ);
  PLA_Finalize( );
  MPI_Finalize( );
  
}
