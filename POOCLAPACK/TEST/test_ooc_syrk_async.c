/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/***************************************************************************
  1 Nov 1999 

  POOCLAPACK: Parallel Out-of-Core Linear Algebra Package

  Out-of-Core extension of PLAPACK: Parallel Linear Algebra Package
		  
  Copyright (c) 1999 Robert van de Geijn and 
  The University of Texas at Austin.

  Unlike PLAPACK, this program is NOT being released as free software; 
  you CANNOT redistribute it and/or modify it without explicit written 
  permission from The University of Texas.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

  Written under the direction of: 
 	Robert van de Geijn, Department of Computer Sciences,
	University of Texas at Austin  78712.    
	rvdg@cs.utexas.edu
 
  POOCLAPACK was written by Wesley Reiley and Robert van de Geijn
 			
***********i****************************************************************/


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
    A_ooc  = NULL, C  = NULL,
    zero   = NULL, one  = NULL;
  int      
    n, k,
    nb_distr, nb_ooc,
    error, parameters, sequential, r12,
    me, nprocs, 
    itype,
    variant, op,
    fd; 
  double 
    time,
    flops;
  char 
    filename[30];

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
    printf("variant? (0 = PAN_PAN, 1 = MAT_PAN):\n");
    scanf("%d", &variant );
    printf("variant = %d\n", variant );
  }
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&parameters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sequential, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&r12, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&variant, 1, MPI_INT, 0, MPI_COMM_WORLD);

  PLA_Set_error_checking( error, parameters, sequential, r12 );

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
      /*
      datatype = MPI_COMPLEX;
      break;
    case 3:
      datatype = MPI_DOUBLE_COMPLEX;
      break;
      */
    default:
      PLA_Abort( "unknown datatype", __LINE__, __FILE__ );
    }

    if (me==0) {
      printf("enter n, k:\n");
      scanf("%d%d", &n, &k );
      printf("n = %d, k = %d\n", n, k );
      printf("enter nb_ooc:\n");
      scanf("%d", &nb_ooc );
      printf("nb_ooc = %d\n", nb_ooc );
    }

    MPI_Bcast(&n,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nb_ooc,1, MPI_INT, 0, MPI_COMM_WORLD);

    sprintf( filename, "file%d", me );
    PLA_File_create( filename, datatype, (int)ceil((n*k)/nprocs), &fd );

    PLA_OOC_Matrix_create_without_fd( datatype, 
		        n, 
		        k,
			templ, 
                        PLA_ALIGN_FIRST, 
                        PLA_ALIGN_FIRST, 
                        &A_ooc ); 

    PLA_Obj_attach_fd( A_ooc, fd, n, TRUE );

    create_problem_ooc( nb_ooc, A_ooc );

    PLA_Matrix_create( datatype, 
		        n,
			n,
		        templ, 
                        PLA_ALIGN_FIRST, 
                        PLA_ALIGN_FIRST, 
                        &C ); 

    PLA_Create_constants_conf_to( A_ooc, NULL, &zero, &one );

    create_problem( C );

    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime ();
      
    PLA_OOC_Syrk_async( nb_ooc, PLA_LOWER_TRIANGULAR, PLA_NO_TRANSPOSE,
		   one, A_ooc, one, C ); 

    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime () - time;

    flops = 1.0 * n * n * k;

    if ( me == 0 ) 
      printf("%d %d time = %f, MFLOPS/node = %10.4lf \n", n, k, time,
	       flops / time * 1.0e-6 / nprocs );

    PLA_File_close( fd );
  }


  PLA_Obj_free( &A_ooc );
  PLA_Obj_free( &C );
  PLA_Obj_free( &zero );
  PLA_Obj_free( &one );

  PLA_Temp_free(&templ);
  PLA_Finalize( );
  MPI_Finalize( );
  
}
