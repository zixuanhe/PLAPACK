/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Spectral_decomp( int uplo, PLA_Obj A, PLA_Obj Q, PLA_Obj diag )
/*
  PLA_Spectral_decomp

  Purpose: Compute spectral decomposition A = Q D Q^T where A is
  a given nxn symmetrix matrix A, D is diagonal with the eigenvalues
  of A on the diagonal and the columns of Q equal the eigenvectors 
  of matrix A.

  input:
  uplo                int
                      indicates whether the symmetric matrix is stored
                      in lower or upper triangular part of A

  A                   PLA_MATRIX
                      matrix to be factored

  Q                   PLA_MATRIX
                      matrix in which to return the eigenvectors

  diag                PLA_MVECTOR of width 1
                      vector in which to return the eigenvalues of A

  output:
  A                   PLA_MATRIX
                      overwritten with junk

  Q                   PLA_MATRIX
                      if Q != NULL columns of Q equal the eigenvectors of A
                      o/w eigenvectors are not computed

  diag                PLA_MVECTOR of width 1
                      eigenvalues of A in ascending order
*/
{
  PLA_Obj
    s = NULL, tridiag = NULL, eigenvalues = NULL, Q_mv = NULL;

  PLA_Template
    templ;
  
  int 
    length;

  MPI_Datatype
    datatype;

  /* Create a vector for the scaling factors computed during the 
     reduction to tridiagonal form and perform reduction */
  PLA_Mvector_create_conf_to( A, 1, &s );

  /* Reduce A to tridiagonal form.  If Q != NULL the Householder transforms
     are accumulated in Q */
  PLA_Tri_red( uplo, A, s, Q );
  
  /* Create a duplicated multiscalar in which to store the
     main diagonal and first subdiagonal of A */
  PLA_Obj_template( A, &templ );
  PLA_Obj_global_length( A, &length );
  PLA_Obj_datatype( A, &datatype );
  PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
		      length, 2, templ, &tridiag );

  /* Copy diagonals to multiscalar */
  PLA_Copy_sym_tridiag_to_msc( uplo, A, tridiag );

  /* Locally solve the tridiagonal eigenproblem */

  /* If Q != NULL, copy Q to a multivector */
  if ( Q != NULL ){
    PLA_Mvector_create_conf_to( Q, length, &Q_mv );
    PLA_Copy( Q, Q_mv );
  }

  PLA_Local_sym_tridiag_eig( tridiag, Q_mv );

  /* Copy the eigenvalues to diag */
  PLA_Obj_vert_split_2( tridiag, 1, &eigenvalues, PLA_DUMMY );
  PLA_Copy( eigenvalues, diag );

  /* if Q != NULL, copy Q_mv back to Q */
  if ( Q != NULL )
    PLA_Copy( Q_mv, Q );


  PLA_Obj_free( &s );
  PLA_Obj_free( &tridiag );
  PLA_Obj_free( &eigenvalues );
  PLA_Obj_free( &Q_mv );
  
  return PLA_SUCCESS;
}

