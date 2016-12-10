/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Apply_sym_House( int uplo, PLA_Obj A, PLA_Obj u, PLA_Obj beta )
/*
  PLA_Apply_sym_House

  Purpose:
  
  Apply a Householder orthogonal similarity transformation to
  matrix A:

  Form    A <- ( I + beta u u^T ) A ( I + beta u u^T )
            =  A + beta u w^T + beta w u^T
   where w = v + beta/2 u^T v u   and  v = A u

   Assumptions:
   A is a PLA_MATRIX, u is a PLA_MVECTOR of width 1 with
   (implicitly) unit first entry, and beta is a PLA_MSCALAR
   duplicated on all nodes.
*/
{
  PLA_Obj
    u_1 = NULL, u_1_copy = NULL,  
    v = NULL, 
    zero = NULL, one = NULL, two = NULL,
    alpha = NULL; 
  double 
    d_two = 2.0;

  PLA_Create_constants_conf_to( A, NULL, &zero, &one );

  /* Set first entry of u to one, saving the old value in u_1_copy */
  PLA_Obj_horz_split_2( u, 1, &u_1,
                               PLA_DUMMY );
  PLA_Mvector_create_conf_to( u_1, 1, &u_1_copy );
  PLA_Local_copy( u_1, u_1_copy );
  PLA_Obj_set_to_one( u_1 );

  /* Compute v = A u */
  PLA_Mvector_create_conf_to( u, 1, &v );
  PLA_Symv( uplo, one, A, u, zero, v );

  /* Compute alpha = u^T v */
  PLA_Mscalar_create_conf_to( beta, PLA_ALL_ROWS, PLA_ALL_COLS, &alpha );
  PLA_Dot( u, v, alpha );

  /* Compute alpha = beta/2 * u^T v */
  PLA_Mscalar_create_conf_to( beta, PLA_ALL_ROWS, PLA_ALL_COLS, &two );
  PLA_Obj_set( two, MPI_DOUBLE, &d_two );

  PLA_Local_scal( beta, alpha );
  PLA_Local_inv_scal( two, alpha );

  /* Compute v = w = v + beta/2 * u^T v u */
  PLA_Local_axpy( alpha, u, v );

  /* Update A = A + beta u w^T + beta w u^T */
  PLA_Syr2( uplo, beta, u, v, A );

  /* Restore first entry of u */
  PLA_Local_copy( u_1_copy, u_1 );

  /* Free temporary objects */
  PLA_Obj_free( &u_1 );
  PLA_Obj_free( &u_1_copy );
  PLA_Obj_free( &v );
  PLA_Obj_free( &zero );
  PLA_Obj_free( &one );
  PLA_Obj_free( &two );
  PLA_Obj_free( &alpha );

  return PLA_SUCCESS;
}
