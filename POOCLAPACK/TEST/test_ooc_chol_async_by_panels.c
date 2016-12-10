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

#if MANUFACTURE == CRAY
#include <mpp/rastream.h>
#endif

#define TILE_NUM 32

int main(int argc, char *argv[])
{
  MPI_Comm 
    comm = MPI_COMM_NULL;
  MPI_Datatype
    datatype;
  PLA_Template 
    templ = NULL;
  PLA_Obj  
    A_ooc[TILE_NUM], C  = NULL, A_in = NULL;
  int      
    n, k,
    nb_distr, nb_ooc, nb_tile, nb_alg,
    error, parameters, sequential, r12,
    me, nprocs, 
    itype,
    variant, op,
    fd[TILE_NUM],
    m_tiles,
    length, local_length, local_width, local_ldim,
    j_tile;

  double 
    time,
    flops;
  char 
    filename[30];

  for(j_tile=0; j_tile < TILE_NUM; j_tile++)
    A_ooc[j_tile] = NULL;

#if MANUFACTURE == CRAY
  set_d_stream( 1 );
#endif

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (me==0) {
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
    printf("turn on r12 checking? (0 = NO, 1 = YES):\n");
    scanf("%d", &r12 );
    printf("r12 checking = %d\n", r12 );
    printf("variant? (0 = PAN_PAN, 1 = MAT_PAN):\n");
    scanf("%d", &variant );
    printf("variant = %d\n", variant );
  }
  MPI_Bcast(&nb_distr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb_alg,   1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&error, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&parameters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sequential, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&r12, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&variant, 1, MPI_INT, 0, MPI_COMM_WORLD);

  pla_Environ_set_nb_alg( PLA_OP_ALL_ALG, nb_alg );

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
      printf("enter n:\n");
      scanf("%d", &n );
      printf("n = %d\n", n );
      printf("enter nb_tile:\n");
      scanf("%d", &nb_tile );
      printf("nb_tile = %d\n", nb_tile );
      printf("enter nb_ooc:\n");
      scanf("%d", &nb_ooc );
      printf("nb_ooc = %d\n", nb_ooc );
    }

    MPI_Bcast(&n,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nb_ooc,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nb_tile,1, MPI_INT, 0, MPI_COMM_WORLD);

    m_tiles = (int)ceil((double)n/(double)nb_tile);


    /* Create array of Matrices */
    for(j_tile=0; j_tile < m_tiles; j_tile++){

      if(j_tile == (m_tiles-1))
	length = n-((m_tiles-1)*nb_tile);
      else
	length = nb_tile;
	  

      PLA_OOC_Matrix_create_without_fd( datatype, 
					length, 
					n,
					templ, 
					j_tile*nb_tile, 
					PLA_ALIGN_FIRST,
					&A_ooc[j_tile] ); 

      sprintf( filename, "file%d_%d", me, j_tile );
      PLA_Obj_local_length( A_ooc[ j_tile ], &local_length );
      PLA_Obj_local_width( A_ooc[ j_tile ], &local_width );
      PLA_Obj_local_ldim(A_ooc[ j_tile ], &local_ldim );
/*      
      PLA_File_create( filename, datatype, (int)ceil((length*n)/nprocs), 
		       &fd[j_tile] );
      PLA_Obj_attach_fd( A_ooc[j_tile], fd[j_tile], length, TRUE );
*/

      PLA_File_create( filename, datatype, local_length*local_width,
		       &fd[j_tile] );
      PLA_Obj_attach_fd( A_ooc[j_tile], fd[j_tile], local_ldim, TRUE );
    }


    PLA_File_stats_start();

    create_problem_ooc_by_panels_spd( nb_tile, m_tiles, A_ooc);

    PLA_File_stats_stop();


    PLA_File_stats_start();

    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime ();
    

    PLA_OOC_Chol_async_by_panels(nb_ooc, nb_tile, m_tiles, A_ooc);


    MPI_Barrier( MPI_COMM_WORLD );
    time = MPI_Wtime () - time;
    PLA_File_stats_stop();

    flops = 1.0/3.0 * n * n * n;

    if ( me == 0 ) {
      printf("%d %d time = %f, MFLOPS/node = %10.4lf \n", n, n, time,
	       flops / time * 1.0e-6 / nprocs );
    }

    for(j_tile=0; j_tile < m_tiles; j_tile++)
      PLA_File_close( fd[j_tile] );
  }
 
  for(j_tile=0; j_tile < TILE_NUM; j_tile++)
    PLA_Obj_free( &A_ooc[j_tile] );
  PLA_Obj_free( &C );

  PLA_Temp_free(&templ);
  PLA_Finalize( );
  MPI_Finalize( );
  
}
