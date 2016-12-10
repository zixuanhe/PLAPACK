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
    m, n,
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
  }
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);

  pla_Environ_set_nb_alg( PLA_OP_ALL_ALG, nb_alg );

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

    if ( me == 0 ){
      printf("enter m, n:\n");
      scanf("%d%d", &m, &n );
      printf("m = %d n = %d\n", m, n );
    }

    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
    PLA_Matrix_create( datatype, 
		       m, 
		       n,
		       templ,
		       PLA_ALIGN_FIRST, 
		       PLA_ALIGN_FIRST, 
		       &A ); 

    create_problem( A );

    PLA_Matrix_create( datatype, 
		       m,
		       n,
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

    PLA_Global_show( "A = [", A, "%6.3lf ", "] " );
    PLA_Global_show( "B = [", B, "%6.3lf ", "] " );

    PLA_Elementwise_mult( A, B, C );

    PLA_Global_show( "C = [", C, "%6.3lf ", "] " );

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
