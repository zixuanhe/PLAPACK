/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Create_svd_problem( int num_transforms, 
			    PLA_Obj diag, PLA_Obj A, PLA_Obj U, PLA_Obj V )
/*
  PLA_Create_svd_problem

  Purpose: Create a matrix with known singular values given in diag.

  input:  num_transformations  indicates the number of random Householder
                               transforms to use to create A.  (An equal
                               number will be applied from left and right).
          diag                 MVECTOR with eigenvalues to be used to create
                               matrix A
  ouput:  A                    symmetrix MATRIX with eigenvalues given in diag
          U                    columns of U equal the left singular vectors of A.
	                       if upon entry U = NULL, U is not computed.
          V                    columns of V equal the right singular vectors of A.
	                       if upon entry V = NULL, V is not computed.

*/
{
  PLA_Obj
    u = NULL, u_1 = NULL, u_1_dup = NULL,
    beta = NULL, temp = NULL, A_diag = NULL;
  int 
    i, proj_old, length_A, width_A;
  double 
    d_minus_two = -2.0;

  /* Set A to diagonal matrix with entries of diag on diagonal */
  PLA_Obj_set_to_zero( A );
  PLA_Obj_global_length( A, &length_A );
  PLA_Obj_global_width ( A, &width_A );
  PLA_Obj_split_4( A, min( length_A, width_A ), min( length_A, width_A ), &A_diag,   PLA_DUMMY, 
		                                                          PLA_DUMMY, PLA_DUMMY );
  PLA_Obj_set_diagonal( A_diag, diag ); 
  
  /* Initialize U = I where U has same dimension as row dimension of A */
  if ( U != NULL )
    PLA_Obj_set_to_identity( U );

  /* Create random Householder transforms and apply to A  from left */

  PLA_Obj_get_orientation( A, &proj_old );
  PLA_Obj_set_orientation( A, PLA_PROJ_ONTO_COL );
  PLA_Mvector_create_conf_to( A, 1, &u );
  PLA_Obj_set_orientation( A, proj_old );
  PLA_Obj_horz_split_2( u, 1, &u_1,
			       PLA_DUMMY );
  PLA_Mscalar_create_conf_to( u_1, PLA_ALL_ROWS, PLA_ALL_COLS, &u_1_dup );
  PLA_Mscalar_create_conf_to( u_1, PLA_ALL_ROWS, PLA_ALL_COLS, &beta );
  PLA_Mscalar_create_conf_to( u_1, PLA_ALL_ROWS, PLA_ALL_COLS, &temp );

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

    /* Apply Householder transform to A from left */
    PLA_Apply_House_transform( PLA_SIDE_LEFT, A, u, beta );

    /* Apply Householder transform to U */
    if ( U != NULL )
      PLA_Apply_House_transform( PLA_SIDE_LEFT, U, u, beta );
  }

  /* Create random Householder transforms and apply to A  from right */

  /* Initialize V = I where V has same dimension as column dimension of A */
  if ( V != NULL )
    PLA_Obj_set_to_identity( V );

  PLA_Obj_get_orientation( A, &proj_old );
  PLA_Obj_set_orientation( A, PLA_PROJ_ONTO_ROW );
  PLA_Mvector_create_conf_to( A, 1, &u );
  PLA_Obj_set_orientation( A, proj_old );
  PLA_Obj_horz_split_2( u, 1, &u_1,
			       PLA_DUMMY );
  PLA_Mscalar_create_conf_to( u_1, PLA_ALL_ROWS, PLA_ALL_COLS, &u_1_dup );
  PLA_Mscalar_create_conf_to( u_1, PLA_ALL_ROWS, PLA_ALL_COLS, &beta );
  PLA_Mscalar_create_conf_to( u_1, PLA_ALL_ROWS, PLA_ALL_COLS, &temp );

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

    /* Apply Householder transform to A from right */
    PLA_Apply_House_transform( PLA_SIDE_RIGHT, A, u, beta );

    /* Apply Householder transform to V */
    if ( V != NULL )
      PLA_Apply_House_transform( PLA_SIDE_LEFT, V, u, beta );
  }

  /* Free temporary objects */
  PLA_Obj_free( &u );
  PLA_Obj_free( &u_1 );
  PLA_Obj_free( &u_1_dup );
  PLA_Obj_free( &beta );
  PLA_Obj_free( &temp );
  PLA_Obj_free( &A_diag );

  return PLA_SUCCESS;
}
    

