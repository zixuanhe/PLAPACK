/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Compute_WY( PLA_Obj A_mv, PLA_Obj s, PLA_Obj W_mv, PLA_Obj Y_mv )

/* to use this routine, we need to modify the calling of
* PLA_QR_right
*/
/*
  Purpose: Compute WY transform from Householder vectors stored
           in A and s.

  Note: Utility routine used as part of computation of QR factorization
        of a matrix.
           
  Input:  A_mv    --   General mxn matrix A   
                       (PLA_MVECTOR)
          s       --   vector of scaling factors
                       (MSCALAR of width=n, duplicated to all nodes)
          W_mv, Y_mv
                  --   matrices for storing W and Y which define the 
                       WY transform.
                       (PLA_MVECTOR of width n)

  Output: W_mv, Y_mv
                  --   W and Y which define the WY transform.

  Assumptions:  n<=m

  Return value: PLA_SUCCESS unless input parameter error is detected.
*/
{
  PLA_Obj   a_B1    = NULL,   A_mv_BR = NULL,   
            s_cur   = NULL,   beta    = NULL,   
            W_L     = NULL,   w_1     = NULL,   W_R = NULL,
            W_BR    = NULL,   w_B1    = NULL,
            Y_BL    = NULL,   y_B1    = NULL,   Y_BR = NULL,
            y_11    = NULL,   u_loc   = NULL,   u    = NULL,
                              u_loc_L = NULL,   u_L  = NULL,
            minus_one = NULL, zero    = NULL,   one     = NULL;
  int       global_width;

  PLA_Create_constants_conf_to( A_mv, &minus_one, &zero, &one );

  PLA_Obj_set_to_zero( W_mv );
  PLA_Obj_set_to_zero( Y_mv );


  /* A_mv_BR tracks the active part of A_mv */
  PLA_Obj_view_all( A_mv, &A_mv_BR );
  /* s_cur   tracks the active part of s */
  PLA_Obj_view_all( s,    &s_cur );

  /* W_L tracks the part of W already computed, W_R the part yet to be 
     computed.  Ditto for W_BR, Y_BL and Y_BR */  
  PLA_Obj_vert_split_2( W_mv, 0,  &W_L,      &W_R );
  PLA_Obj_vert_split_2( W_mv, 0,  PLA_DUMMY, &W_BR );
  PLA_Obj_vert_split_2( Y_mv, 0,  &Y_BL,     &Y_BR );

  /* Create duplicated multiscalar to hold u and local contributions to u */
  PLA_Mscalar_create_conf_to( s, PLA_ALL_ROWS, PLA_ALL_COLS, &u );
  PLA_Mscalar_create_conf_to( s, PLA_ALL_ROWS, PLA_ALL_COLS, &u_loc );

  /* u_L and u_loc_L tracks the part of u and u_loc 
     corresponding to Y_BL and Y_BR, respectively */
  PLA_Obj_horz_split_2( u, 0,      &u_L,     PLA_DUMMY );
  PLA_Obj_horz_split_2( u_loc, 0,  &u_loc_L, PLA_DUMMY );

  while( TRUE ){
    PLA_Obj_global_width( A_mv_BR, &global_width );
    if ( 0 == global_width ) break;

    /* Split off next column of A_mv_BR and columns of W and Y 
       to be computed */
    PLA_Obj_vert_split_2( A_mv_BR,  1,     &a_B1, PLA_DUMMY );
    PLA_Obj_vert_split_2( W_R,      1,     &w_1,  &W_R );
    PLA_Obj_vert_split_2( W_BR,     1,     &w_B1, PLA_DUMMY );
    PLA_Obj_vert_split_2( Y_BR,     1,     &y_B1, PLA_DUMMY );

    /* Split off current scaling factor beta */
    PLA_Obj_horz_split_2( s_cur,    1,     &beta,
                                           &s_cur );

    /* y_B1 = a_B1 with first element set to 1 */
    PLA_Local_copy( a_B1, y_B1 );
    PLA_Obj_horz_split_2( y_B1, 1, &y_11, 
                                   PLA_DUMMY );
    PLA_Obj_set_to_one( y_11 );

    /* w_1 = - beta ( /   0   \
                    \  y_B1 / + W_L Y_BL y_B1 ) */
    PLA_Local_copy( y_B1, w_B1 );

    PLA_Local_gemv( PLA_TRANS, one, Y_BL, y_B1, zero, u_loc_L );

    PLA_Reduce( u_loc_L, MPI_SUM, u_L );

    PLA_Obj_global_width( W_L, &global_width );
    if ( global_width > 0 )
      PLA_Local_gemv( PLA_NO_TRANS, one, W_L, u_L, one, w_1 );

    PLA_Local_scal( beta, w_1 );
/*    PLA_Local_scal( minus_one, w_1 ); */

    /* Update views */
    PLA_Obj_view_shift( W_L,    0,
                            0,      1,
                                0 );

    PLA_Obj_view_shift( Y_BL,    1,
                            0,        1,
                                 0 );


    PLA_Obj_view_shift( u_L,     0,
                            0,        0,
                                 1 );

    PLA_Obj_view_shift( u_loc_L,     0,
                               0,        0,
                                     1 );
    
    PLA_Obj_split_4( A_mv_BR, 1, 1,   PLA_DUMMY, PLA_DUMMY,
                                      PLA_DUMMY, &A_mv_BR);
    PLA_Obj_split_4( Y_BR, 1, 1,      PLA_DUMMY, PLA_DUMMY,
                                      PLA_DUMMY, &Y_BR );
    PLA_Obj_split_4( W_BR,    1, 1,   PLA_DUMMY, PLA_DUMMY,
                                      PLA_DUMMY, &W_BR );
  }

  /* Clean up temporary objects */
  PLA_Obj_free( &a_B1 );    PLA_Obj_free( &A_mv_BR );   
  PLA_Obj_free( &s_cur );  PLA_Obj_free( &beta );   
  PLA_Obj_free( &W_L );    PLA_Obj_free( &w_1 );   PLA_Obj_free( &W_R );
  PLA_Obj_free( &W_BR );    PLA_Obj_free( &w_B1 );
  PLA_Obj_free( &Y_BL );   PLA_Obj_free( &y_B1 );   PLA_Obj_free( &Y_BR );
  PLA_Obj_free( &y_11 ); 
  PLA_Obj_free( &u );    PLA_Obj_free( &u_loc ); 
  PLA_Obj_free( &u_L );    PLA_Obj_free( &u_loc_L ); 
  PLA_Obj_free( &minus_one );
  PLA_Obj_free( &zero );   PLA_Obj_free( &one );

  return PLA_SUCCESS;
}  

    
