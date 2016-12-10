/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_SVD( PLA_Obj A, PLA_Obj U, PLA_Obj D, PLA_Obj V )
/*
  PLA_SVD

  Purpose: Compute singular value decomposition A = U D V^T where A is
  a given matrix A, D is diagonal with the singular values
  of A on the diagonal and the columns of U and V equal the left- and
  right-singular vectors, respectively.

  input:
  A                   PLA_MATRIX
                      matrix to be factored

  U                   PLA_MATRIX
                      matrix in which to return the left singular vectors

  D                   PLA_MVECTOR of width 1
                      vector in which to return the singular values of A

  V                   PLA_MATRIX
                      matrix in which to return the right singular vectors

  output:
  A                   PLA_MATRIX
                      overwritten with junk

  U                   PLA_MATRIX
                      if U != NULL columns of V equal the left singular vectors of A
                      o/w not computed

  diag                PLA_MVECTOR of width 1
                      eigenvalues of A in ascending order

  V                   PLA_MATRIX
                      if V != NULL columns of V equal the right singular vectors of A
                      o/w not computed
*/
{
  PLA_Obj
    sL = NULL, sR = NULL, bidiag = NULL, singularvalues = NULL, U_mv = NULL, V_mv = NULL;

  PLA_Template
    templ;
  
  int 
    length_A, width_A, value;

  MPI_Datatype
    datatype;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_SVD_enter( A, U, D, V );

  /* Create vectors for the scaling factors computed during the 
     reduction to bidiagonal form and perform reduction */
  PLA_Mvector_create_conf_to( A, 1, &sL );

  PLA_Obj_set_orientation( A, PLA_PROJ_ONTO_ROW );
  PLA_Mvector_create_conf_to( A, 1, &sR );
  PLA_Obj_set_orientation( A, PLA_PROJ_ONTO_COL );

  /* Reduce A to bidiagonal form.  If U != NULL and/or V != NULL the Householder transforms
     are accumulated in U and V */
  PLA_Bi_red( PLA_UPPER_TRIANGULAR, A, sL, sR, U, V );

  /* Create a duplicated multiscalar in which to store the
     main diagonal and first superdiagonal of A */
  PLA_Obj_template( A, &templ );
  PLA_Obj_global_length( A, &length_A );
  PLA_Obj_global_width( A, &width_A );
  PLA_Obj_datatype( A, &datatype );
  PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
		      min( length_A, width_A ), 2, templ, &bidiag );

  /* Copy diagonals to multiscalar */
  PLA_Copy_bidiag_to_msc( PLA_UPPER_TRIANGULAR, A, bidiag );

  /* Locally solve the tridiagonal eigenproblem */

  /* If U != NULL, copy U to a multivector */
  if ( U != NULL ){
    PLA_Mvector_create_conf_to( U, length_A, &U_mv );
    PLA_Copy( U, U_mv );
  }

  /* If V != NULL, copy V^T to a multivector */
  if ( V != NULL ){
    PLA_Mvector_create_conf_to( V, width_A, &V_mv );
    PLA_Copy( V, V_mv );
  }

  PLA_Local_bidiag_svd( bidiag, U_mv, V_mv  );

  /* Copy the eigenvalues to D */
  PLA_Obj_vert_split_2( bidiag, 1, &singularvalues, PLA_DUMMY );
  PLA_Copy( singularvalues, D );

  /* if U != NULL, copy U_mv back to U */
  if ( U != NULL )
    PLA_Copy( U_mv, U );

  /* if V != NULL, copy V_mv back to V */
  if ( V != NULL )
    PLA_Copy( V_mv, V );


  PLA_Obj_free( &sL );
  PLA_Obj_free( &sR );
  PLA_Obj_free( &bidiag );
  PLA_Obj_free( &singularvalues );
  PLA_Obj_free( &U_mv );
  PLA_Obj_free( &V_mv );
  

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_SVD_exit( A, U, D, V );

  return PLA_SUCCESS;
}

