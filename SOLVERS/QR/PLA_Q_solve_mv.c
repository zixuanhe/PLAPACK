/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Q_solve_mv( int side, int trans, PLA_Obj A, PLA_Obj s, PLA_Obj b)
/*
  Purpose: Solve one of the following:

       OPERATION             SIDE               TRANS
       ----------------------------------------------
       Q   X = B               PLA_SIDE_LEFT      PLA_NO_TRANS
       Q^T X = B               PLA_SIDE_LEFT      PLA_TRANS
       X Q   = B               PLA_SIDE_RIGHT     PLA_NO_TRANS
       X Q^T = B               PLA_SIDE_RIGHT     PLA_TRANS
       
    where Q is stored as Householder vectors below the diagonal of A and in s
    and B is given.  The solution X overwrites B.

  Input:  side    --   Indicates whether Q appears to left or right of X
                       (int)
          trans   --   Indicates whether to transpose Q
                       (int)
          A       --   General mxn matrix A   
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)
	  b       --   right-hand-side(s)
                       (MVECTOR)

  Output: b       --   first n elements contain least-square solution

  Return value: PLA_SUCCESS iff QR factorization is completed successfully
*/
{
  if ( side == PLA_SIDE_LEFT ) {
    if ( trans == PLA_NO_TRANS )
      return PLA_Q_solve_mv_l_n( A, s, b );
    else if ( trans == PLA_TRANS ) {
/*      return PLA_Q_solve_mv_l_t( A, s, b ); */
      printf("PLA_Q_solve_mv: left, trans not yet implemented\n");
      exit( 0 );
    }
    else {
      printf("PLA_Q_solve_mv: illegal value for trans\n");
      exit( 0 );
    }
  }
  else if ( side == PLA_SIDE_RIGHT ) {
    if ( trans == PLA_NO_TRANS ) {
/*      return PLA_Q_solve_mv_r_n( A, s, b ); */
      printf("PLA_Q_solve_mv: right, notrans not yet implemented\n");
      exit( 0 );
    }
    else if ( trans == PLA_TRANS ) {
/*      return PLA_Q_solve_mv_r_t( A, s, b ); */
      printf("PLA_Q_solve_mv: right, trans not yet implemented\n");
      exit( 0 );
    }
    else {
      printf("PLA_Q_solve_mv: illegal value for trans\n");
      exit( 0 );
    }
  }
  else {
    printf("PLA_Q_solve_mv: illegal value for side\n");
    exit( 0 );
  }
}


int PLA_Q_solve_mv_l_n( PLA_Obj A, PLA_Obj s, PLA_Obj b)
/*
  Purpose: Solve Q x = b where Q is stored as Householder vectors
           below the diagonal of A and in s.

  Input:  A       --   General mxn matrix A   
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)
	  b       --   right-hand-side(s)
                       (MVECTOR)

  Output: A       --   QR factorization.  R is stored in upper-triangular
                       portion of A.  Q is stored in vector form below
                       the diagonal of A, with scaling factors in vector s
          s       --   vector of scalars for Householder transforms
	  b       --   first n elements contain least-square solution

  Return value: PLA_SUCCESS iff QR factorization is completed successfully
*/
{
  int       size, me, size_in, dummy;
  PLA_Obj   ABR         = NULL,
            A1          = NULL,     A2              = NULL,
            A1_mv       = NULL,     A1_mvBR         = NULL,     
	    u           = NULL,     u1             = NULL,
            sB          = NULL,     s1             = NULL,
            s1_dup      = NULL,     s1_dup_B       = NULL,
            beta        = NULL,
            bB          = NULL,     alpha          = NULL;
  PLA_Template  templ;
  MPI_Datatype datatype;

  PLA_Obj_template( A, &templ );
  PLA_Obj_datatype( A, &datatype );
  PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
		      1, 1, templ, &alpha );

  /* Initially partition A = / ATL |  *  \  where
                             \ ABL | ABR /  ATL is 0x0  */
  PLA_Obj_view_all( A, &ABR );

  /* Initially partition s = / sT \
                             \ sB / where sB is of length 0 */
  PLA_Obj_view_all( s, &sB );

  /* Initially partition b = / bT \
                             \ bB / where bT is of length 0 */
  PLA_Obj_view_all( b, &bB );

  while ( TRUE ) {
    /* Determine block size */
    PLA_Obj_split_size( ABR, PLA_SIDE_LEFT, &size, &dummy );
    if ( size == 0 ) break;

    /* Partition ABR = ( A1, A2 ) */
    PLA_Obj_vert_split_2( ABR, size,    &A1, &A2 );

    /* Partition sB  = ( s1, 
                         sB ) */
    PLA_Obj_horz_split_2( sB, size, &s1, 
                                    &sB );
    /* Create a duplicated copy of s1, to duplicate  beta's */
    PLA_Mscalar_create_conf_to( s1, PLA_ALL_ROWS, PLA_ALL_COLS, 
                                &s1_dup );

    /* Copy the current panel to a multivector */
    PLA_Mvector_create_conf_to( A1, size, &A1_mv );

    PLA_Copy( s1, s1_dup );
    PLA_Copy( A1, A1_mv );

    PLA_Obj_view_all( A1_mv, &A1_mvBR );
    PLA_Obj_view_all( s1_dup, &s1_dup_B );

    while ( TRUE ){
      PLA_Obj_split_size( A1_mvBR, PLA_SIDE_LEFT, &size_in, &dummy );
      if ( size_in == 0 ) break;

      PLA_Obj_vert_split_2( A1_mvBR, 1, &u, PLA_DUMMY );
      PLA_Obj_horz_split_2( u, 1,   &u1, 
                                    PLA_DUMMY );
      PLA_Obj_set_to_one( u1 );

      PLA_Obj_horz_split_2( s1_dup_B, 1, &beta, 
                                         &s1_dup_B );
      PLA_Dot( u, bB, alpha );
      PLA_Local_scal( beta, alpha );
      PLA_Local_axpy( alpha, u, bB );
      
      PLA_Obj_horz_split_2( bB, 1, PLA_DUMMY,
			           &bB );
      PLA_Obj_split_4( A1_mvBR, 1, 1, PLA_DUMMY, PLA_DUMMY,
                                      PLA_DUMMY, &A1_mvBR );
    }

    /* Update view ABR to view currently active part of A */
    PLA_Obj_horz_split_2( A2, size,          PLA_DUMMY,
                                             &ABR );
  }

  /* Free temporary objects and views */
  PLA_Obj_free( &ABR ); 
  PLA_Obj_free( &A1 );            PLA_Obj_free( &A2 );
  PLA_Obj_free( &A1_mv );         PLA_Obj_free( &A1_mvBR );         
  PLA_Obj_free( &u );             PLA_Obj_free( &u1 );
  PLA_Obj_free( &sB );            PLA_Obj_free( &s1 );
  PLA_Obj_free( &s1_dup );        PLA_Obj_free( &s1_dup_B );
  PLA_Obj_free( &beta );
  PLA_Obj_free( &bB );            PLA_Obj_free( &alpha );

  return( PLA_SUCCESS );
}


