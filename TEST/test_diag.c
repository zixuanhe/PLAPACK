/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

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
    A      = NULL, x    = NULL, B = NULL;
  int      
    m, n,
    nb_distr, 
    error, parameters, sequential, 
    me, nprocs, 
    itype,
    iside, side;
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
  }
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&parameters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sequential, 1, MPI_INT, 0, MPI_COMM_WORLD);

  PLA_Set_error_checking( error, parameters, sequential, FALSE );

  PLA_Comm_1D_to_2D_ratio(MPI_COMM_WORLD, 1.0, &comm);
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
      printf("enter m, n:\n");
      scanf("%d %d", &m, &n );
      printf("m = %d n = %d\n", m, n );
      printf("enter side (0=left, 1=right):\n");
      scanf("%d", &iside );
      printf("iside = %d\n", iside );
    }

    MPI_Bcast(&m,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iside, 1, MPI_INT, 0, MPI_COMM_WORLD);

    switch( iside ){
    case 0: 
      side = PLA_SIDE_LEFT;
      break;
    case 1: 
      side = PLA_SIDE_RIGHT;
      break;
    default:
      PLA_Abort( "side case not supported", __LINE__, __FILE__ );
    }

    PLA_Matrix_create( datatype, 
		        n, n,
			templ, 
                        PLA_ALIGN_FIRST, 
                        PLA_ALIGN_FIRST, 
                        &A ); 

    create_problem( A );

    PLA_Mvector_create( datatype, 
		        n,
			1,
			templ, 
			PLA_ALIGN_FIRST, 
			&x ); 
      
    create_problem( x );

    {
      PLA_Obj
	A_BR = NULL,   A_11 = NULL;
      int
	size;

      double i = 1.0;

      PLA_Obj_view_all( A, &A_BR );

      while ( TRUE ){
	PLA_Obj_global_length( A_BR, &size );
	if ( size == 0 ) break;

	PLA_Obj_split_4( A_BR, 1, 1, &A_11,     PLA_DUMMY,
			             PLA_DUMMY, &A_BR );

	PLA_Obj_set( A_11, MPI_DOUBLE, (void * ) &i );
	
	i = i+1.0;
      }
      
      PLA_Obj_free( &A_BR );
      PLA_Obj_free( &A_11 );
    }
	  
    PLA_Global_show( "A = ", A, "%5.2lf ", " " );
    PLA_Global_show( "x before= ", x, "%5.2lf ", " " );
    PLA_Multiply_by_diagonal( side, A, x );
    PLA_Global_show( "x after = ", x, "%5.2lf ", " " );

    
    PLA_Matrix_create( datatype, 
		       n, m,
		       templ, 
		       PLA_ALIGN_FIRST, 
		       PLA_ALIGN_FIRST, 
		       &B ); 

    create_problem( B );

    PLA_Global_show( "A = ", A, "%5.2lf ", " " );
    PLA_Global_show( "B before= ", B, "%5.2lf ", " " );
    PLA_Multiply_by_diagonal( PLA_SIDE_LEFT, A, B );
    PLA_Global_show( "B after = ", B, "%5.2lf ", " " );

    PLA_Matrix_create( datatype, 
		       m, n,
		       templ, 
		       PLA_ALIGN_FIRST, 
		       PLA_ALIGN_FIRST, 
		       &B ); 

    create_problem( B );

    PLA_Global_show( "A = ", A, "%5.2lf ", " " );
    PLA_Global_show( "B before= ", B, "%5.2lf ", " " );
    PLA_Multiply_by_diagonal( PLA_SIDE_RIGHT, A, B );
    PLA_Global_show( "B after = ", B, "%5.2lf ", " " );
  }

  PLA_Obj_free( &A );
  PLA_Obj_free( &B );
  PLA_Obj_free( &x );

  PLA_Temp_free(&templ);
  PLA_Finalize( );
  MPI_Finalize( );
  
}
