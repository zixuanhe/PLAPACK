/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Tri_red( int uplo, PLA_Obj A, PLA_Obj s, PLA_Obj Q )
/*
  PLA_Tri_red

  Purpose: Reduce symmetric matrix A to tridiagonal form using 
  Householder similarity transformations.

  input:
  uplo                    indicates whether A is stored in upper or
                          lower triangular part 
  A                       MATRIX to be reduced

  output:
  A                       Reduced matrix A.  Householder vectors used
                          to reduce A are stored below first subdiagonal
                          of A.
  s                       Scaling factors for the Householder transforms
                          computed to reduce A.
  Q                       if Q != NULL, Q equals the accumulation of 
                          Householder transforms.
*/
{
  PLA_Obj
    u = NULL,  u_B = NULL, 
    beta_B = NULL, beta_1 = NULL, beta_1_dup = NULL,
    A_BR = NULL, a_21 = NULL, A_21 = NULL, 
    q_11 = NULL, q_12 = NULL, q_21 = NULL, Q_22 = NULL;
  int
    size, value = PLA_SUCCESS;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Tri_red_enter( uplo, A, s, Q );

  if ( uplo != PLA_LOWER_TRIANGULAR ) 
    PLA_Abort( "only uplo == PLA_LOWER_TRIANGULAR currently supported",
		__LINE__, __FILE__ );

  /* Create a vector in which to compute the Householder vector */
  PLA_Mvector_create_conf_to( A, 1, &u );
  
  /* Create a duplicated multiscalar in which to hold the scaling factor
     for the Householder transform being computed */
  PLA_Obj_horz_split_2( s, 1, &beta_1,
			      PLA_DUMMY );
  PLA_Mscalar_create_conf_to( beta_1, PLA_ALL_ROWS, PLA_ALL_COLS, 
			      &beta_1_dup );

  /* Track the active parts of A, s, and u */
  PLA_Obj_view_all( A, &A_BR );    
  PLA_Obj_view_all( s, &beta_B );    
  PLA_Obj_view_all( u, &u_B );

  while ( TRUE ){  
    PLA_Obj_global_length( A_BR, &size );
    if ( 1 == size ) break;

    /* Partition A_BR = / alpha_11    *    \
                        \   a_21    A_BR   /  where alpha_11 is 1x1 */
    PLA_Obj_split_4( A_BR, 1, 1, PLA_DUMMY, PLA_DUMMY,
    		                 &a_21,     &A_BR );
    /* Split of the current element of vector s */
    PLA_Obj_horz_split_2( beta_B, 1, &beta_1,
                                     &beta_B );
    /* View the part of u in which to compute the Householder vector */
    PLA_Obj_horz_split_2( u_B, 1, PLA_DUMMY,
			          &u_B );

    /* Redistributed a_21 as a vector and compute Householder transform */
    PLA_Copy( a_21, u_B );  
    PLA_Compute_House_v( u_B, beta_1_dup );

    /* Place data back in A and s */
    PLA_Copy( u_B, a_21 );
    PLA_Local_copy( beta_1_dup, beta_1 );

    /* Update A_BR <- ( I - beta_1 u_B u_B^T ) A_BR ( I - beta_1 u_B u_B^T ) */
    PLA_Apply_sym_House( uplo, A_BR, u_B, beta_1_dup );
  }

  if ( Q != NULL ){
    /* Compute the orthogonal matrix */
    PLA_Obj_split_4( Q, 1, 1, &q_11, &q_12,
		              &q_21, &Q_22 );

    PLA_Obj_set_to_one ( q_11 );
    PLA_Obj_set_to_zero( q_12 );
    PLA_Obj_set_to_zero( q_21 );

    PLA_Obj_split_4( A, 1, -1, PLA_DUMMY, PLA_DUMMY,
		                &A_21,     PLA_DUMMY );

    PLA_Form_Q( PLA_NO_TRANSPOSE, A_21, s, Q_22 );
  }

  /* Free the temporary objects */
  PLA_Obj_free( &u );
  PLA_Obj_free( &u_B );
  PLA_Obj_free( &beta_B );
  PLA_Obj_free( &beta_1 );
  PLA_Obj_free( &beta_1_dup );
  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &a_21 );
  PLA_Obj_free( &A_21 );
  PLA_Obj_free( &q_11 );
  PLA_Obj_free( &q_12 );
  PLA_Obj_free( &q_21 );
  PLA_Obj_free( &Q_22 );

  if ( PLA_ERROR_CHECKING )   
    value = PLA_Tri_red_exit( uplo, A, s, Q );

  return PLA_SUCCESS;
}
