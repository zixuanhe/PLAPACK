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
    quit,
    length_top, width_left, width_right, length_bottom,
    nb_distr, 
    error, parameters, sequential, r12,
    me, nprocs, 
    itype;

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
    printf("turn on error checking? (0 = NO, 1 = YES):\n");
    scanf("%d", &error );
    printf("error checking = %d\n", error );
    printf("turn on parameter checking? (0 = NO, 1 = YES):\n");
    scanf("%d", &parameters );
    printf("parameter checking = %d\n", parameters );
    printf("turn on sequential checking? (0 = NO, 1 = YES):\n");
    scanf("%d", &sequential );
    printf("sequential checking = %d\n", sequential );
    printf("turn on r12 checking? (0 = NO, 1 = YES):\n");
    scanf("%d", &r12 );
    printf("r12 checking = %d\n", r12 );
  }
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&parameters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sequential, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&r12, 1, MPI_INT, 0, MPI_COMM_WORLD);

  PLA_Set_error_checking( error, parameters, sequential, r12 );

  PLA_Comm_1D_to_2D_ratio(MPI_COMM_WORLD, 1.0, &comm);
  PLA_Init( comm );
    
  PLA_Temp_create( nb_distr, 0, &templ );
    
  if (me==0) {
    printf("enter object type:\n");
    printf(" 0 = matrix\n");
    printf(" 1 = mvector\n");
    printf(" 2 = mscalar\n");
    printf(" 3 = pmvector\n");
    scanf("%d", &itype );
    printf("itype = %d\n", itype );
    printf("enter m, n:\n");
    scanf("%d%d", &m, &n );
  }
  MPI_Bcast(&itype, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  switch( itype ){
  case 0:
    PLA_Matrix_create( MPI_DOUBLE, m, n, templ, 
                        PLA_ALIGN_FIRST, 
                        PLA_ALIGN_FIRST, 
                        &A ); 
    break;
  case 1:
    PLA_Mvector_create( 
			MPI_DOUBLE, m, n, templ, 
                        PLA_ALIGN_FIRST, 
                        &A ); 
    break;
  case 2:
    PLA_Mscalar_create( 
			MPI_DOUBLE, 
			PLA_ALL_ROWS, PLA_ALL_COLS,
			m, n, templ, 
                        &A ); 
    break;
  case 3:
    PLA_Pmvector_create( 
			 MPI_DOUBLE, PLA_PROJ_ONTO_COL,
			 PLA_ALL_COLS, m, n, templ, 
			 PLA_ALIGN_FIRST, 
			 &A ); 
    break;
  default:
    PLA_Abort( "unknown datatype", __LINE__, __FILE__ );
  }

  while ( TRUE ){
    if (me==0) {
      printf("quit? ( 1 = yes 0 = no )\n");
      scanf("%d", &quit );
    }
    MPI_Bcast(&quit, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if ( quit == 1 ) break;

    if (me==0) {
      printf("enter length_top, etc.:\n");
      scanf("%d%d%d%d", &length_top, &width_left, &width_right,
	                &length_bottom );
      printf("length top = %d\n", length_top );
      printf("width left  = %d\n", width_left );
      printf("width right = %d\n", width_right );
      printf("length top = %d\n", length_bottom );
    }
    MPI_Bcast(&length_top, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&width_left, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&width_right, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&length_bottom, 1, MPI_INT, 0, MPI_COMM_WORLD);

    PLA_Obj_view_shift( A,   length_top,
		       width_left, width_right,
		              length_bottom );

  }
  PLA_Obj_free(&A );

  PLA_Temp_free(&templ);
  PLA_Finalize( );
  MPI_Finalize( );
  
}
