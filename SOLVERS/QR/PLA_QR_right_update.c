/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

g#include "PLA.h"

int PLA_QR_right_update( PLA_Obj A, PLA_Obj s )

/*
  Purpose: Compute Householder transform based QR factorization of A

  Input:  A       --   General mxn matrix A   
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)

  Output: A       --   QR factorization.  R is stored in upper-triangular
                       portion of A.  Q is stored in vector form below
                       the diagonal of A, with scaling factors in vector s
          s       --   vector of scalars for Householder transforms

  Return value: PLA_SUCCESS iff QR factorization is completed successfully
*/
{
  int       size, me, width_A,
            nb_alg1, nb_alg2, nb_alg;
  PLA_Template  templ;
  PLA_Obj   ABR         = NULL,
            A1          = NULL,     A2              = NULL,
            A1_mv       = NULL,     
            W_mv        = NULL,     
            Y_mv        = NULL,     
            sB          = NULL,     s1              = NULL,
            s1_dup      = NULL;

  PLA_Obj_global_width( A, &width_A );

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
    /* Create a duplicated copy of s1, in which to compute the
       beta's */
    PLA_Mscalar_create_conf_to( s1, PLA_ALL_ROWS, PLA_ALL_COLS, 
                                &s1_dup );

    /* Copy the current panel to a multivector */
    PLA_Mvector_create_conf_to( A1, size, &A1_mv );

    PLA_Copy( A1, A1_mv );

    /* Factor the current panel */
    {
      PLA_Obj 
	A1_mv_T = NULL,
	A1_mv_B = NULL;
      int width, length;

      PLA_Obj_global_width( A1_mv, &width );
      PLA_Obj_horz_split_2( A1_mv, width, &A1_mv_T,
			                  &A1_mv_B );

      PLA_Obj_horz_split_2( A1_mv, -width_A, PLA_DUMMY,
			                     &A1_mv_B );

      PLA_QR_mv_update( A1_mv_T, A1_mv_B, s1_dup );
      /* PLA_QR_mv( A1_mv, s1_dup ); */
    }

    /* Place factored panel back in the right place and place
       the computed beta's in correct part of s */
    PLA_Copy( A1_mv, A1 );
    PLA_Copy( s1_dup, s1 );

    /* Create multivectors for W and Y and compute WY transform */
    PLA_Mvector_create_conf_to( A1, size, &W_mv );
    PLA_Mvector_create_conf_to( A1, size, &Y_mv );
    /*    PLA_Compute_WY_update( A1_mv, s1_dup, W_mv, Y_mv ); */
    {
      PLA_Obj 
	A1_mv_T = NULL,
	A1_mv_B = NULL,
	W1_mv_B = NULL,
	Y1_mv_B = NULL;
      int width, length;

      PLA_Obj_global_width( A1_mv, &width );
      PLA_Obj_horz_split_2( A1_mv, width, &A1_mv_T,
			                 PLA_DUMMY );

      PLA_Obj_horz_split_2( A1_mv, -width_A, PLA_DUMMY,
			                     &A1_mv_B );

      PLA_Obj_horz_split_2( W_mv, -width_A, PLA_DUMMY,
			                     &W1_mv_B );

      PLA_Obj_horz_split_2( Y_mv, -width_A, PLA_DUMMY,
			                     &Y1_mv_B );

      PLA_Compute_WY_update2( A1_mv_T, A1_mv_B, s1_dup, W1_mv_B, Y1_mv_B );
      /* PLA_QR_mv( A1_mv, s1_dup ); */
    }

    /* A2 <- ( I + W Y^T ) A2 */
    {
      PLA_Obj
	A_T = NULL, A_B = NULL, Y_B = NULL, W_B = NULL;
      int width;
      
      PLA_Obj_global_width( A1_mv, &width );

      PLA_Obj_horz_split_2( A2, width, &A_T,
			              PLA_DUMMY );

      PLA_Obj_horz_split_2( A2, -width_A, PLA_DUMMY,
			                  &A_B );

      PLA_Obj_horz_split_2( W_mv, -width_A, PLA_DUMMY,
			                    &W_B );

      PLA_Obj_horz_split_2( Y_mv, -width_A, PLA_DUMMY,
			                    &Y_B );

      pla_Apply_W_Y_transform_L_T_update( W_B, Y_B, A_T, A_B );
    
      PLA_Obj_free( &A_T );
      PLA_Obj_free( &A_B );
      PLA_Obj_free( &W_B );
      PLA_Obj_free( &Y_B );
    }

    /* Update view ABR to view currently active part of A */
    PLA_Obj_horz_split_2( A2, size,          PLA_DUMMY,
                                             &ABR );
  }

  /* Free temporary objects and views */
  PLA_Obj_free( &ABR ); 
  PLA_Obj_free( &A1 );            PLA_Obj_free( &A2 );
  PLA_Obj_free( &A1_mv );         
  PLA_Obj_free( &W_mv );          PLA_Obj_free( &Y_mv );          
  PLA_Obj_free( &sB );            PLA_Obj_free( &s1 );
  PLA_Obj_free( &s1_dup );

  return( PLA_SUCCESS );
}


