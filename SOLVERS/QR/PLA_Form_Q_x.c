/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Form_Q_x( int trans, int uplo, PLA_Obj A, PLA_Obj s, PLA_Obj Q )

/*
  Purpose: Compute Q = H_1 * H_2 * ... * H_n-1 where H_i is
  the ith Householder transform stored in the ith column of A of
  uplo == PLA_LOWER_TRIANGULAR or the ith row of A if 
  uplo == PLA_UPPER_TRIANGULAR and 
  ith entry in s.
  If trans = PLA_NO_TRANS, Q is computed.  Otherwise, Q^T is computed.

  Input:  trans   --   Indicates whether Q or Q^T is to be computed.
          uplo    --   Indicates whether Householder vectors are
                       stored above or below diagonal.
          A       --   General mxn matrix A 
                       (PLA_MATRIX)
          s       --   vector for storing scalar in Householder transforms
                       (MVECTOR of length=min(m,n), width=1)

  Output: Q      

  Return value: PLA_SUCCESS iff Q is computed successfully.
*/
{
  int       size, length, width, nb_alg1, nb_alg2, nb_alg, proj_old;
  PLA_Template  templ;
  PLA_Obj   ATL         = NULL,     ABR            = NULL,
            QBR         = NULL,
            A1B         = NULL,     AB1_mv         = NULL,     
            W_mv        = NULL,     Y_mv        = NULL,     
            sL          = NULL,
            s_cur       = NULL,     s_dup   = NULL;

  if ( uplo == PLA_LOWER_TRIANGULAR )
    return PLA_Form_Q( trans, A, s, Q );

  /* Determine algorithmic blocking size */
  PLA_Obj_template( A, &templ );
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg1 ); 
  PLA_Environ_nb_alg( PLA_OP_PAN_MAT, templ, &nb_alg2 ); 
  nb_alg = ( nb_alg1 > nb_alg2 ? nb_alg1 : nb_alg2 );

  /* Q = I */
  PLA_Obj_set_to_identity( Q );

  /* Apply the Householder vectors like 
                 ( H_1 ( ... ( H_n-1 ( H_n I ) ) ... ) ) */

  /* Initially partition A = / ATL | ATR \  where ATL is square and
                             \  *  | ABR /  ABR has width 0  */

  PLA_Obj_global_length( A, &length );
  PLA_Obj_global_width ( A, &width );
  width = min( length, width );
  PLA_Obj_split_4( A, width, width,   &ATL,      PLA_DUMMY,
                                      PLA_DUMMY, &ABR );

  PLA_Obj_split_4( Q, width, width,   PLA_DUMMY, PLA_DUMMY,
                                      PLA_DUMMY, &QBR );

  PLA_Obj_horz_split_2( s, width,     &sL,       
                                      PLA_DUMMY );

  while ( TRUE ) {
    /* Check if done */

    PLA_Obj_global_length( ATL, &length ); 
    PLA_Obj_global_width( ATL, &width ); 
    if ( ( size = min( min( length, width), nb_alg ) ) == 0 ) break;

    /* Grow ABR by the panel from which to compute the next 
       WY transform.  Grow QBR similarly */

    PLA_Obj_view_shift( ABR,    -size, 
                     -size,               0,
                                  0 );

    PLA_Obj_view_shift( QBR,      -size,
                          -size,           0,
                                    0 );

    /* Partition off A1B from which to compute the next WY transform */
    PLA_Obj_horz_split_2( ABR, size, &A1B, 
                                     PLA_DUMMY );
    PLA_Obj_set_orientation( A1B, PLA_PROJ_ONTO_ROW );

    /* Partition off the scaling factors and duplicate to all nodes */

    PLA_Obj_horz_split_2( sL,  -size,     &sL, 
                                          &s_cur );

    PLA_Mscalar_create_conf_to( s_cur, PLA_ALL_ROWS, PLA_ALL_COLS, 
                                &s_dup );

    /* Redistribute A1B as a multivector and compute W and Y */

    PLA_Obj_get_orientation( QBR, &proj_old );
    PLA_Obj_set_orientation( QBR, PLA_PROJ_ONTO_ROW );
    PLA_Mvector_create_conf_to( QBR, size, &AB1_mv );
    PLA_Mvector_create_conf_to( QBR, size, &W_mv );
    PLA_Mvector_create_conf_to( QBR, size, &Y_mv );
    PLA_Obj_set_orientation( QBR, proj_old );

    PLA_Copy( A1B, AB1_mv );
    PLA_Copy( s_cur, s_dup );

    PLA_Compute_WY( AB1_mv, s_dup, W_mv, Y_mv );

    if ( trans == PLA_NO_TRANSPOSE ) {
      /* Update QBR <- ( I + W Y^T ) QBR */
      PLA_Apply_W_Y_transform ( PLA_SIDE_LEFT, PLA_NO_TRANSPOSE,
                                W_mv, Y_mv, QBR ); 
    }
    else {
      /* Update QBR <- QBR ( I + W^T Y ) QBR */
      PLA_Apply_W_Y_transform ( PLA_SIDE_RIGHT, PLA_TRANSPOSE,
                                W_mv, Y_mv, QBR ); 
    }

    /* Update view of ATL */

    PLA_Obj_split_4( ATL, -size, -size,  &ATL,      PLA_DUMMY,
                                         PLA_DUMMY, PLA_DUMMY );
  }

  PLA_Obj_free( &ATL );           PLA_Obj_free( &ABR );            
  PLA_Obj_free( &QBR ); 
  PLA_Obj_free( &A1B );           PLA_Obj_free( &AB1_mv );
  PLA_Obj_free( &W_mv );          PLA_Obj_free( &Y_mv );          
  PLA_Obj_free( &s_cur );         PLA_Obj_free( &s_dup );
  PLA_Obj_free( &sL );

  return PLA_SUCCESS;
}
