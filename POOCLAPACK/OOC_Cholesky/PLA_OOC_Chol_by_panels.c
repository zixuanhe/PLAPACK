/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_OOC_Chol_by_panels( int nb_ooc, int nb_tile, int m_tiles, 
			     PLA_Obj *A_row_panels )
{

  PLA_Obj
    L_10 = NULL, A_11 = NULL, A_11_in = NULL, 
    L_20_1 = NULL, A_21_1 = NULL, A_21_1_in = NULL,
    min_one = NULL, one = NULL, temp=NULL;

  int 
    size, tile_size,  value,
    j_tile, i_tile;

/*  if ( PLA_ERROR_CHECKING )
    value = PLA_OOC_Chol_enter( nb_ooc, nb_tile, A ); */
  
  for ( j_tile=0; j_tile<m_tiles; j_tile++ ){

    PLA_Create_constants_conf_to(A_row_panels[ j_tile ], 
				  &min_one, NULL, &one);

    PLA_Obj_global_length(A_row_panels[ j_tile ], &tile_size);

    size = j_tile*nb_tile;

    PLA_Obj_vert_split_2( A_row_panels[ j_tile ], size, 
			   &L_10, &temp );
    
    PLA_Obj_vert_split_2( temp, tile_size, &A_11, PLA_DUMMY);

    PLA_Matrix_create_conf_to(A_11, &A_11_in);
 
    PLA_Copy(A_11, A_11_in);

    PLA_OOC_Syrk(nb_ooc, PLA_LOWER_TRIANGULAR, PLA_NO_TRANS,
			min_one, L_10, one, A_11_in);

    PLA_Chol(PLA_LOWER_TRIANGULAR, A_11_in);

    PLA_Copy(A_11_in, A_11);
    
    for ( i_tile=j_tile+1; i_tile<m_tiles; i_tile++ ){
      PLA_Obj_vert_split_2( A_row_panels[ i_tile ], size,
			     &L_20_1, &temp );
      PLA_Obj_vert_split_2( temp, tile_size, &A_21_1, PLA_DUMMY );

      PLA_Matrix_create_conf_to(A_21_1, &A_21_1_in);
      PLA_Copy(A_21_1, A_21_1_in);

      PLA_OOC_Gemm(nb_ooc, PLA_NO_TRANS, PLA_TRANS, min_one, L_20_1, 
		    L_10, one, A_21_1_in);

      PLA_Trsm(PLA_SIDE_RIGHT, PLA_LOWER_TRIANGULAR, PLA_TRANS,
                PLA_NONUNIT_DIAG, one, A_11_in, A_21_1_in);


      PLA_Copy(A_21_1_in, A_21_1);
    }
  }

  PLA_Obj_free( &L_10 );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_11_in );
  PLA_Obj_free( &L_20_1 );
  PLA_Obj_free( &A_21_1 );
  PLA_Obj_free( &A_21_1_in );
  PLA_Obj_free( &min_one );
  PLA_Obj_free( &one );
  PLA_Obj_free( &temp );


/*  if ( PLA_ERROR_CHECKING )   
    value = PLA_OOC_Chol_exit( nb_ooc, nb_tile, A ); */

}

