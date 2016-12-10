/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Bi_red( PLA_Obj A, PLA_Obj sL, PLA_Obj sR, PLA_Obj U, PLA_Obj V )
/*
  PLA_Tri_red

  Purpose: Reduce general matrix A to bidiagonal form using 
  Householder similarity transformations.

  input:
  A                       MATRIX to be reduced

  output:
  A                       Reduced matrix A.  Householder vectors used
                          to reduce A are stored below subdiagonal
                          of A and above first superdiagonal of A.
  sL                      Scaling factors for the Householder transforms
                          computed to reduce A (applied to left of A).
  sR                      Scaling factors for the Householder transforms
                          computed to reduce A (applied to right of A).
  U                       if U != NULL, U equals the accumulation of 
                          Householder transforms that were applied from the
                          left.
  V                       if V != NULL, V equals the accumulation of 
                          Householder transforms that were applied from the
                          right.
*/
{
  PLA_Obj
    uv = NULL,  uv_B = NULL, u = NULL,  v = NULL,
    sR_B = NULL, betaR_1 = NULL, betaR_1_dup = NULL,
    sL_B = NULL, betaL_1 = NULL, betaL_1_dup = NULL,
    A_BR = NULL, a_21 = NULL, A_21 = NULL, 
    u_11 = NULL, u_12 = NULL, u_21 = NULL, U_22 = NULL;
    v_11 = NULL, v_12 = NULL, v_21 = NULL, V_22 = NULL;
  int
    size, value = PLA_SUCCESS;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Bi_red_enter( A, sR, sL, U, V );

  /* The whole purpose of the game is to update most of the matrix
     A = A - ( u | w_R ) / 0 | w_L^T \
                         \ 0 | v^T   /
     where u equals the Householder transform that hits
     the active part of A from the left and
     v equals the Householder transform that hits the active
           part from the right:
     A <- ( I - betaL_1 u u^T ) A ( I - betaR_1 / 0 \ ( 0 v^T ) 
                                                \ v /                 */

  /* Create a multivector in which to compute ( u  |  w_R ) */

  PLA_Mvector_create_conf_to( A, 1, &uwR );

  /* Create a multivector in which to compute ( w_L | v ) */

  PLA_Mvector_create_conf_to( A, 1, &wLv );
  
  /* Create duplicated multiscalars in which to hold the scaling factor
     for the Householder transforms being computed */
  PLA_Obj_horz_split_2( sL, 1, &betaL_1,
			      PLA_DUMMY );
  PLA_Mscalar_create_conf_to( betaL_1, PLA_ALL_ROWS, PLA_ALL_COLS, 
			      &betaL_1_dup );
  PLA_Obj_horz_split_2( sR, 1, &betaR_1,
			      PLA_DUMMY );
  PLA_Mscalar_create_conf_to( betaR_1, PLA_ALL_ROWS, PLA_ALL_COLS, 
			      &betaR_1_dup );

  /* Track the active parts of A, sL, sR, and uwR and wLv */
  PLA_Obj_view_all( A, &A_BR );    
  PLA_Obj_view_all( sL, &betaL_B );    
  PLA_Obj_view_all( sR, &betaR_B );    
  PLA_Obj_view_all( uwR, &uwR_B );
  PLA_Obj_horz_split_2( wLv, 1, PLA_DUMMY,
                                &wLv_B );

  while ( TRUE ){  
    PLA_Obj_global_length( A_BR, &size );
    if ( 0 == size ) break;

    /* Partition A_BR = ( a_L      A_R    ) where a_L is a column from
       which to compute the next u */
    PLA_Obj_vert_split_2( A_BR, 1, &a_L, &A_R );
    /* Split of the current element of vector sL */
    PLA_Obj_horz_split_2( betaL_B, 1, &betaL_1,
                                      &betaL_B );
    /* View the part of uwR in which to compute u and wR */
    PLA_Obj_vert_split_2( uwR_B, 1,  &u, &wR );

    /* View the part of wLv in which to compute w_L and v */
    PLA_Obj_vert_split_2( wLv_B, 1,  &w_L, &v );

    /* Redistributed a_L as a vector and compute Householder transform */
    PLA_Copy( a_L, u );  
    PLA_Compute_House_v( u, betaL_1_dup );

    /* Place data back in A and sL */
    PLA_Copy( u, a_L );
    PLA_Local_copy( betaL_1_dup, betaL_1 );

    /* Partition A_BR = /  *    a_12^T \  
                        \  *    A_22   / splitting off first column and row */
    PLA_Obj_vert_split_4( A_BR, 1, 1, PLA_DUMMY, &a_12,
                                      PLA_DUMMY, &A_BR );
    
    /* Copy a_12 into v */
    PLA_Obj_set_orientation( PLA_PROJ_ONTO_ROW, a_12 );
    PLA_Copy( a_12, v );
    
    /* w_L = betaL_1 * A_BR^T * u */
       First must split off first element of u */
    PLA_Obj_horz_split_2( u, 1, PLA_DUMMY,
			        &u_2 );

    PLA_Gemv( PLA_TRANSPOSE, betaL_1_dup, A_BR, u_2, zero, w_L );

    /* In v update a_12 = a_12 - betaL_1 * A_R^T u = 
                        = a_12 - betaL_1 * a_12 - w_L */

    PLA_Negate( betaL_1 );
    PLA_Axpy( betaL_1, v, v );
    PLA_Axpy( minus_one, w_L, v );

    /* Compute Householder transform to reflect resulting v */
    PLA_Compute_House_v( v, betaR_1_dup );

    /* Place data back in A and sR */
    PLA_Copy( v, a_12 );
    PLA_Local_copy( betaR_1_dup, betaR_1 );

    /* Set first element of v = 1 */
    PLA_Obj_horz_split_2( v, 1, &v_1,
			        PLA_DUMMY );
    PLA_Obj_set_to_one( v_1 );

    /* Update A_BR <- ( I - betaL_1 u_2 u_2^T ) A_BR ( I - betaR_1 v v^T )
                   =  ( A_BR - u_2 ( betaL_1 u_2^T A_BR ) ) ( I - betaR_1 v v^T )
		   =  A_BR - ( u_2 | w_R ) * ( w_L | v )
       where w_L = betaL_1 * A_BR^T u_2,  w_R = betaR_1 * A_BR * v - betaR_1 * w_L^T * v * u_2 
    */
    
    /* Compute w_R */
    PLA_Gemv( PLA_NO_TRANSPOSE, betaR_1, A_BR, v, zero, w_R );
    PLA_Dot( x_L, v, temp );
    PLA_Scal( betaR_1_dup, temp );
    PLA_Negate( temp );
    PLA_Axpy( temp, u_2, w_R );

    /* Update A_BR */
    PLA_Pmvector_create_conf_to( A_BR, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 2, &uwR_2_dup );
    PLA_Copy( uwR_2, uwR_2_dup );

    PLA_Pmvector_create_conf_to( A_BR, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 2, &wLv_dup );
    PLA_Copy( wLv, wLv_dup );
    
    PLA_Gemm( PLA_NO_TRANSPOSE, PLA_NO_TRANSPOSE, minus_one, uwR_2, wLv, one, A_BR );
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
    value = PLA_Bi_red_exit( A, sL, sR, U, V );

  return PLA_SUCCESS;
}
