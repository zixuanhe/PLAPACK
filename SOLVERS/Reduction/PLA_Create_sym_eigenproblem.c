/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Create_sym_eigenproblem( int uplo, int num_transforms, 
				  PLA_Obj diag, PLA_Obj A, PLA_Obj Q )
/*
  PLA_Create_sym_eigenproblem

  Purpose: Create a symmetric matrix, filling only the part indicated
  by uplo, with known eigenvalues given in diag.

  input:  uplo                 indicates whether to fill upper or lower
                               triangular part of A
          num_transformations  indicates the number of random Householder
                               similarity transforms to us to create A
          diag                 MVECTOR with eigenvalues to be used to create
                               matrix A
  ouput:  A                    symmetrix MATRIX with eigenvalues given in diag
          Q                    columns of A equal the eigenvectors of A
	                       if upon entry Q = NULL, eigenvectors are not
                               computed
*/
{
  PLA_Obj
    u = NULL, u_1 = NULL, u_1_dup = NULL,
    beta = NULL, temp = NULL;
  int 
    i;
  double 
    d_minus_two = -2.0;

  /* Set A to diagonal matrix with entries of diag on diagonal */
  PLA_Obj_set_to_zero( A );
  PLA_Obj_set_diagonal( A, diag );
  
  /* Initialize Q = I */
  if ( Q != NULL )
    PLA_Obj_set_to_identity( Q );

  PLA_Mvector_create_conf_to( A, 1, &u );
  PLA_Obj_horz_split_2( u, 1, &u_1,
			       PLA_DUMMY );
  PLA_Mscalar_create_conf_to( u_1, PLA_ALL_ROWS, PLA_ALL_COLS, &u_1_dup );
  PLA_Mscalar_create_conf_to( u_1, PLA_ALL_ROWS, PLA_ALL_COLS, &beta );
  PLA_Mscalar_create_conf_to( u_1, PLA_ALL_ROWS, PLA_ALL_COLS, &temp );

  /* Create random Householder transforms and apply to A */
  for ( i=0; i<num_transforms; i++ ){
    /* Create a random vector with entries ( -1,1 ) */
    create_random_data( u );

    /* Normalize the vector to have unit first entry */
    PLA_Copy( u_1, u_1_dup );
    PLA_Local_inv_scal( u_1_dup, u );

    /* Compute beta = - 2.0/ u^t u */
    PLA_Obj_set( beta, MPI_DOUBLE, &d_minus_two );
    PLA_Dot( u, u, temp );
    PLA_Local_inv_scal( temp, beta );

    /* Apply Householder similarity transform to A */
    PLA_Apply_sym_House( uplo, A, u, beta );

    /* Apply Householder transform to Q */
    if ( Q != NULL )
      PLA_Apply_House_transform( PLA_SIDE_LEFT, Q, u, beta );
  }

  /* Free temporary objects */
  PLA_Obj_free( &u );
  PLA_Obj_free( &u_1 );
  PLA_Obj_free( &u_1_dup );
  PLA_Obj_free( &beta );
  PLA_Obj_free( &temp );

  return PLA_SUCCESS;
}
    

