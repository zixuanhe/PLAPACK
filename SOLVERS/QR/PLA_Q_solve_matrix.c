/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Q_solve_matrix( PLA_Obj A, PLA_Obj s, PLA_Obj B )
/*
  Purpose: Solve Q X = B where Q is stored as Householder vectors
           below the diagonal of A and in s.

  Input:  A       --   General mxn matrix A   
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)
	  B       --   right-hand-side(s)
                       (MVECTOR or MATRIX of length m)

  Output: A       --   QR factorization.  R is stored in upper-triangular
                       portion of A.  Q is stored in vector form below
                       the diagonal of A, with scaling factors in vector s
          s       --   vector of scalars for Householder transforms
	  B       --   first n rows contain least-square solution

  Return value: PLA_SUCCESS iff QR factorization is completed successfully
*/
{
  int       size, me,
            nb_alg1, nb_alg2, nb_alg;
  PLA_Template  templ;
  PLA_Obj   ABR         = NULL,
            A1          = NULL,     A2              = NULL,
            A1_mv       = NULL,     
            W_mv        = NULL,     
            Y_mv        = NULL,     
            sB          = NULL,     s1              = NULL,
            s1_dup      = NULL,     BB              = NULL;

  /* Determine algorithmic blocking size */
  PLA_Obj_template( A, &templ );
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg1 ); 
  PLA_Environ_nb_alg( PLA_OP_PAN_MAT, templ, &nb_alg2 ); 
  nb_alg = ( nb_alg1 > nb_alg2 ? nb_alg1 : nb_alg2 );

  /* Initially partition A = / ATL |  *  \  where
                             \ ABL | ABR /  ATL is 0x0  */
  PLA_Obj_view_all( A, &ABR );

  /* Initially partition s = / sT \
                             \ sB / where sB is of length 0 */
  PLA_Obj_view_all( s, &sB );

  /* Initially partition B = / BT \  where
                             \ BB /  BT has 0 rows */
  PLA_Obj_view_all( B, &BB );

  while ( TRUE ) {
    /* Check if done */
    PLA_Obj_global_width( ABR , &size ); 
    if ( ( size = min( size, nb_alg ) ) == 0 ) break;

    /* Partition ABR = ( A1, A2 ) */
    PLA_Obj_vert_split_2( ABR,      size, &A1, &A2 );

    /* Partition sB  = ( s1, 
                         sB ) */
    PLA_Obj_horz_split_2( sB, size, &s1, 
                                    &sB );


    /* Create a duplicated copy of s1, in which to duplicate
       beta's */
    PLA_Mscalar_create_conf_to( s1, PLA_ALL_ROWS, PLA_ALL_COLS, 
                                &s1_dup );
    PLA_Copy( s1, s1_dup );

    /* Copy the current panel to a multivector */
    PLA_Mvector_create_conf_to( A1, size, &A1_mv );

    PLA_Copy( A1, A1_mv );

    /* Create multivectors for W and Y and compute WY transform */
    PLA_Mvector_create_conf_to( A1, size, &W_mv );
    PLA_Mvector_create_conf_to( A1, size, &Y_mv );
    PLA_Compute_WY( A1_mv, s1_dup, W_mv, Y_mv );

    /* BB <- ( I + W Y^T ) BB */
    PLA_Apply_W_Y_transform ( PLA_SIDE_LEFT, PLA_TRANSPOSE,
                              W_mv, Y_mv, BB); 

    /* Update view ABR to view currently active part of A */
    PLA_Obj_horz_split_2( A2, size,          PLA_DUMMY,
                                             &ABR );

    /* Update view BB to view currently active part of B */
    PLA_Obj_horz_split_2( BB, size,          PLA_DUMMY,
                                             &BB );
  }

  /* Free temporary objects and views */
  PLA_Obj_free( &ABR ); 
  PLA_Obj_free( &A1 );            PLA_Obj_free( &A2 );
  PLA_Obj_free( &A1_mv );         
  PLA_Obj_free( &W_mv );          PLA_Obj_free( &Y_mv );          
  PLA_Obj_free( &sB );            PLA_Obj_free( &s1 );
  PLA_Obj_free( &s1_dup );
  PLA_Obj_free( &BB );

  return( PLA_SUCCESS );
}


