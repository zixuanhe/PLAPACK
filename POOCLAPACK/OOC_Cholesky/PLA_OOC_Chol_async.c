/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"


int PLA_OOC_Chol_async(int nb_ooc, int nb_tile, PLA_Obj A )
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

  PLA_Obj_split_4(A, 0, 0, PLA_DUMMY, PLA_DUMMY,
  	                    &L_BL,      &A_BR);

  PLA_Create_constants_conf_to(A, &min_one, NULL, &one);

  while(TRUE) {
    PLA_Obj_global_length(A_BR, &size);
    if (size == 0) break;
    size = min(nb_tile, size);

    PLA_Obj_horz_split_2(L_BL, size, &L_10,
		                      &L_20);

    PLA_Obj_split_4(A_BR, size, size, &A_11, PLA_DUMMY,
		                       &A_21, &A_BR);

    PLA_Matrix_create_conf_to(A_11, &A_11_in);
    PLA_Copy(A_11, A_11_in);

    PLA_OOC_Syrk_async(nb_ooc, PLA_LOWER_TRIANGULAR, PLA_NO_TRANS,
		  min_one, L_10, one, A_11_in);

    PLA_Chol(PLA_LOWER_TRIANGULAR, A_11_in);

    PLA_Copy(A_11_in, A_11);
    
    PLA_Obj_view_all(L_20, &L_20_B);
    PLA_Obj_view_all(A_21, &A_21_B);

    while (TRUE) {
      PLA_Obj_global_length(L_20_B, &size_in);
      if(size_in == 0) break;
      size_in = min(nb_tile, size_in);

      PLA_Obj_horz_split_2(L_20_B, size_in, &L_20_1,
		                             &L_20_B);

      PLA_Obj_horz_split_2(A_21_B, size_in, &A_21_1,
		                             &A_21_B);

      PLA_Matrix_create_conf_to(A_21_1, &A_21_1_in);
      PLA_Copy(A_21_1, A_21_1_in);

      PLA_OOC_Gemm_async(nb_ooc, PLA_NO_TRANS, PLA_TRANS, min_one, L_20_1, 
		    L_10, one, A_21_1_in);
      
      PLA_Trsm(PLA_SIDE_RIGHT, PLA_LOWER_TRIANGULAR, PLA_TRANS,
                PLA_NONUNIT_DIAG, one, A_11_in, A_21_1_in);


      PLA_Copy(A_21_1_in, A_21_1);
    }


    PLA_Obj_view_shift(L_BL,     size, 
                               0,       size, 
                                   0);


  }

  PLA_Obj_free( &L_BL );
  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &L_10 );
  PLA_Obj_free( &L_20 );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_21 );
  PLA_Obj_free( &A_11_in );
  PLA_Obj_free( &L_20_B );
  PLA_Obj_free( &A_21_B );
  PLA_Obj_free( &L_20_1 );
  PLA_Obj_free( &A_21_1 );
  PLA_Obj_free( &A_21_1_in );
  PLA_Obj_free( &min_one );
  PLA_Obj_free( &one );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_OOC_Chol_exit( nb_ooc, nb_tile, A );

}

