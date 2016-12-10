/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_OOC_Chol(int nb_ooc, int nb_tile, PLA_Obj A )
/*
  Purpose : Parallel Out-of-Core Cholesky Factorization.  
            This particular version assumes only the lower triangular 
            portion of A is stored.

  IN     nb_ooc      integer, out-of-core algorithmic block size
                     (typically setting this to the in-core algorithmic
                      block size for Cholesky factorization works great.
                      E.g., set nb_ooc = 128)	    
  IN     nb_tile     integer, tile size.  This should be set to a size that
                     allows two tiles to reside in the combined memories of
                     the processors, leaving some room for a few smaller
                     nb_tile x nb_ooc panels.
  IN/OUT A           matrix to be factored, store out on disk.

  Uses left-looking variant of Cholesky factorization:

  ******************************************************************
      
      Partition  A = / L_TL ||   *  \
                     | =====  ===== |
                     \ L_BL || A_BR /
     	     where L_TL is 0 x 0 
      while A_BR is not 0 x 0 
         Determine block size t
	 Partition  A = / L_TL ||   *  \    / L_00 ||  *   |  *   \ 
                        | =====  ===== | =  | ====    ====   ==== |
                        \ L_BL || A_BR /    | L_10 || A_11 |  *   |
                                            | ----    ----   ---- |
                                            \ L_20 || A_21 | A_22 /
     	        where A_00 = A_TL and A_11 is t x t 
         Update A_11 <- A_11 - L_10 * L_10' (out-of-core syrk update)
         Update A_11 <- L_11 = Chol. Fact.( A_11 )
         Update A_21 <- A_21 - L_20 * L_10' (out-of-core gemm)
	 Update A_21 <- L_21 =  A_21 * inv(L_11') (out-of-core trsm)
         Continue with
	            A = / L_TL ||   *  \    / L_00 |  *   ||  *   \ 
                        | =====  ===== | =  | ----   ----    ---- |
                        \ L_BL || A_BR /    | L_10 | L_11 ||  *   |
                                            | ====   ====    ==== |
                                            \ L_20 | L_21 || A_22 /
      endwhile
              
  ******************************************************************

  NOTE:  For details on how to implement OOC parallel Cholesky factorization
         see
	 
	 Wesley Reiley and Robert van de Geijn,
	 "POOCLAPACK: Parallel Out-of-Core Linear Algebra Package,"
	 submitted to ICS2000.

	 Wesley Reiley, 
         "..."
	 undergraduate honor's thesis.
*/
{

  PLA_Obj
    L_BL = NULL, A_BR = NULL, L_10 = NULL, L_20 = NULL,
    A_11 = NULL, A_21 = NULL, A_11_in = NULL, 
    L_20_B = NULL, A_21_B = NULL, L_20_1 = NULL,
    A_21_1 = NULL, A_21_1_in = NULL,
    min_one = NULL, one = NULL;

  int 
    size, size_in,  value, me, nprocs, width, length;

  double
    flops, time;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_OOC_Chol_enter( nb_ooc, nb_tile, A );

  /* Create usual duplicated scalar constants */
  PLA_Create_constants_conf_to(A, &min_one, NULL, &one);

  /*
      Partition  A = / L_TL ||   *  \
                     | =====  ===== |
                     \ L_BL || A_BR /
		     where A_TL is 0 x 0          
      Note: L_TL is not needed for the computation, so we don't view it  */

  PLA_Obj_split_4(A, 0, 0, PLA_DUMMY, PLA_DUMMY,
  	                    &L_BL,      &A_BR);


  while(TRUE) {
    /* Check how big to pick the next A_11 */
    PLA_Obj_global_length(A_BR, &size);
    if (size == 0) break;
    size = min(nb_tile, size);

    /*
	 Partition  A = / L_TL ||   *  \    / L_00 ||  *   |  *   \ 
                        | =====  ===== | =  | ====    ====   ==== |
                        \ L_BL || A_BR /    | L_10 || A_11 |  *   |
                                            | ----    ----   ---- |
                                            \ L_20 || A_21 | A_22 / */
    PLA_Obj_horz_split_2(L_BL, size, &L_10,
		                      &L_20);

    PLA_Obj_split_4(A_BR, size, size, &A_11, PLA_DUMMY,
		                       &A_21, &A_BR);

    /* Create an incore matrix into which to copy OOC matrix A_11 */
    PLA_Matrix_create_conf_to(A_11, &A_11_in);

    /* Bring A_11 in-core */
    PLA_Copy(A_11, A_11_in);

    /* Update in-core A_11 <- A_11 - L_10 * L_10' */
    PLA_OOC_Syrk(nb_ooc, PLA_LOWER_TRIANGULAR, PLA_NO_TRANS,
		  min_one, L_10, one, A_11_in);

    /* Factor in-core A_11 */
    PLA_Chol(PLA_LOWER_TRIANGULAR, A_11_in);

    /* Write updated in-core A_11 out to disk */
    PLA_Copy(A_11_in, A_11);
    
    /* Update A_21 <- A_21 - L_20 * L_10' followed by
              A_21 <- L_21 = A_21 * inv(L_11') 

       Algorithm:
 
       *****************************************************************

       Partition  L_20 = / L_20_T \    A_21 = / A_21_T \
                         | ====== |           | ====== |
                         \ L_20_B /           \ A_21_B /
     	     where L_20_T has 0 rows and A_21_T has 0 rows.
       while L_20_B has rows left
	 Partition  / L_20_T \   / L_20_0 \     / A_21_T \   / A_21_0 \
                    | ====== | = | ====== | and | ====== | = | ====== |    
                    \ L_20_B /   | L_20_1 |     \ A_21_B /   | A_21_1 |
                                 | ------ |                  | ------ |
                                 \ L_20_2 /                  \ A_21_2 /
    	        where L_20_0 = L_20_T, L_20_1 has length t, and
                      A_21_0 = A_21_T, A_21_1 has length t    
         Update A_21_1 <- A_21_1 - L_20_1 * L_10' (out-of-core gemm)
	 Update A_21_1 <- L_21_1 = A_21_1 * inv(L_11') (in-core trsm)
         Continue with
	            / L_20_T \   / L_20_0 \     / A_21_T \   / A_21_0 \
                    | ====== | = | ------ | and | ====== | = | ------ |
                    \ L_20_B /   | L_20_1 |     \ A_21_B /   | A_21_1 |
                                 | ====== |                  | ====== |    
                                 \ L_20_2 /                  \ A_21_2 /
       endwhile

       *****************************************************************
    */              

    PLA_Obj_view_all(L_20, &L_20_B);
    PLA_Obj_view_all(A_21, &A_21_B);

    while (TRUE) {
      /* Determine next height of L_20_1 and A_21_1 */
      PLA_Obj_global_length(L_20_B, &size_in);
      if(size_in == 0) break;
      size_in = min(nb_tile, size_in);

      /*
	 Partition  / L_20_T \   / L_20_0 \     / A_21_T \   / A_21_0 \
                    | ====== | = | ====== | and | ====== | = | ====== |    
                    \ L_20_B /   | L_20_1 |     \ A_21_B /   | A_21_1 |
                                 | ------ |                  | ------ |
                                 \ L_20_2 /                  \ A_21_2 /
      */
      PLA_Obj_horz_split_2(L_20_B, size_in, &L_20_1,
		                             &L_20_B);

      PLA_Obj_horz_split_2(A_21_B, size_in, &A_21_1,
		                             &A_21_B);

      /* Bring A_21_1 in-core */
      PLA_Matrix_create_conf_to(A_21_1, &A_21_1_in);
      PLA_Copy(A_21_1, A_21_1_in);

      /* Update A_21_1 <- A_21_1 - L_20_1 * L_10' (out-of-core gemm) */
      PLA_OOC_Gemm(nb_ooc, PLA_NO_TRANS, PLA_TRANS, min_one, L_20_1, 
		    L_10, one, A_21_1_in);
      
      /* Update A_21_1 <- L_21_1 = A_21_1 * inv(L_11') (in-core trsm) */
      PLA_Trsm(PLA_SIDE_RIGHT, PLA_LOWER_TRIANGULAR, PLA_TRANS,
                PLA_NONUNIT_DIAG, one, A_11_in, A_21_1_in);

      /* Write A_21_1 back to disk */
      PLA_Copy(A_21_1_in, A_21_1);
    }

    /* Update view of L_BL */
    PLA_Obj_view_shift(L_BL,     size, 
                               0,       size, 
                                   0);
  }

  /* free temporary objects and views */
  PLA_Obj_free( &L_BL );         PLA_Obj_free( &A_BR );
  PLA_Obj_free( &L_10 );         PLA_Obj_free( &L_20 );
  PLA_Obj_free( &A_11 );         PLA_Obj_free( &A_21 );
  PLA_Obj_free( &A_11_in );      PLA_Obj_free( &L_20_B );
  PLA_Obj_free( &A_21_B );       PLA_Obj_free( &L_20_1 );
  PLA_Obj_free( &A_21_1 );       PLA_Obj_free( &A_21_1_in );
  PLA_Obj_free( &min_one );      PLA_Obj_free( &one );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_OOC_Chol_exit( nb_ooc, nb_tile, A );

  return value;
}

