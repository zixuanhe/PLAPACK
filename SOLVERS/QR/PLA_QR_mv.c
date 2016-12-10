/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_QR_mv( PLA_Obj A_mv, PLA_Obj s )

/*
  Purpose: Compute Householder transform based QR factorization of A

  Note: Utility routine used as part of computation of QR factorization
        of a matrix.
           
  Input:  A_mv    --   General mxn matrix A   
                       (PLA_MVECTOR)
          s       --   vector for storing scalar in Householder transforms
                       (MSCALAR of width=n, duplicated to all nodes)

  Output: A_mv    --   QR factorization.  R is stored in upper-triangular
                       portion of A.  Q is stored in vector form below
                       the diagonal of A, with scaling factors in vector s
          s       --   vector of scalars for Householder transforms

  Assumptions:  n <= m

  Return value: PLA_SUCCESS iff QR factorization is completed successfully
*/
{
  int   size;
  PLA_Template  templ;
  PLA_Obj A_mv_cur = NULL,      v_1 = NULL,
          a1 = NULL,            v = NULL,
          s_cur = NULL,         beta = NULL, 
	  u = NULL,             u_loc = NULL,
	  zero = NULL,          one = NULL;
  MPI_Datatype datatype;

  PLA_Obj_template( A_mv, &templ );
  PLA_Obj_datatype( A_mv, &datatype );

  /* Create usual duplicated scalar constants */
  PLA_Create_constants_conf_to( A_mv, NULL, &zero, &one);

  /* View all of A_mv and s */
  PLA_Obj_view_all(A_mv,&A_mv_cur);
  PLA_Obj_view_all(s,&s_cur);

  /* Create a vector for copy of current column of A_mv */
  PLA_Mvector_create_conf_to( A_mv, 1, &v );

  /* create 2 multiscalar to hold result of local part of 
     trans(A_mv_cur) * v and global result of same operation */
  PLA_Obj_global_width(A_mv, &size);
  PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
                      size, 1, templ, &u_loc );
  PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
                      size, 1, templ, &u );
  while ( TRUE) {
    PLA_Obj_global_width( A_mv_cur, &size );
    if ( 0 == size ) break;
    /* Split off first column of A_mv_cur */
    PLA_Obj_vert_split_2( A_mv_cur, 1,     &a1, &A_mv_cur);
    /* Split off first element of s_cur to store scaling factor beta */
    PLA_Obj_horz_split_2( s_cur, 1,        &beta,
			            	   &s_cur);
    /* Compute v and beta so that ( I - beta v v^T ) a1 = +- || a1 ||_2  
       where first element of v is unity, and rest of v overwrites 
       corresponding part of a1 */

    PLA_Compute_House_v( a1, beta); 

    if ( 1 == size ) break;

    /* let v = a1, with first element replaced by 1 */
    PLA_Local_copy( a1, v );
    PLA_Obj_horz_split_2(v, 1,     &v_1,
			           PLA_DUMMY);
    PLA_Obj_set_to_one( v_1 );

    /* size u and u_loc */
    PLA_Obj_horz_split_2(u,     1,    PLA_DUMMY,
			              &u );
    PLA_Obj_horz_split_2(u_loc, 1,    PLA_DUMMY,
			              &u_loc );
    /* compute locally u_loc  = beta*trans(A)*a1 */
    PLA_Local_gemv( PLA_TRANSPOSE, beta, A_mv_cur, v, zero, u_loc ); 
    /* reduce to all into the duplicated multiscalar u */
    PLA_Reduce( u_loc, MPI_SUM, u );

    /* update A_mv_cur = A_mv_cur + v u^T */
    PLA_Local_ger(one, v, u, A_mv_cur);

    /* update the splittings */
    PLA_Obj_horz_split_2( A_mv_cur, 1, PLA_DUMMY,
			               &A_mv_cur);
    PLA_Obj_horz_split_2( v, 1,        PLA_DUMMY,
			               &v );
  }
  /* some clean up */
  PLA_Obj_free( &A_mv_cur );       PLA_Obj_free( &v_1 );
  PLA_Obj_free( &a1 );             PLA_Obj_free( &v );
  PLA_Obj_free( &s_cur );          PLA_Obj_free( &beta );
  PLA_Obj_free( &u );              PLA_Obj_free( &u_loc );
  PLA_Obj_free( &zero );           PLA_Obj_free( &one );

  return( PLA_SUCCESS );
}

