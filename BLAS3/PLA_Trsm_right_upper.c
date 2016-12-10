/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Trsm_right_upper( int transa, int diag,
			  PLA_Obj alpha, PLA_Obj A, PLA_Obj B )
{
  PLA_Obj 
    A_22 = NULL,                     A_22_dup = NULL,
    A_12 = NULL, A_TL = NULL,        A_12t_dup = NULL,  
    B_L = NULL,  B_2 = NULL,         C_2_dup = NULL,
    one = NULL,  minus_one = NULL;

  int 
    nb_alg, size;

  PLA_Template
    templ;

  PLA_Local_scal( alpha, B );   
  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );
  
  PLA_Obj_global_width( B, &size );
  nb_alg = min( size, nb_alg );

  if ( transa == PLA_NO_TRANS ){
    PLA_Abort( "case not yet implemented", __LINE__, __FILE__ );
  }
  else { /* transa == PLA_TRANS */
    PLA_Obj_template( A, &templ );
    PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

    PLA_Obj_view_all( A, &A_TL );        PLA_Obj_view_all( B, &B_L );

    PLA_Pmvector_create_conf_to( 
  	  B, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, nb_alg, &C_2_dup );
    PLA_Pmvector_create_conf_to( 
          B, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, nb_alg, &A_12t_dup );

    while ( TRUE ){
      PLA_Obj_global_length( A_TL, &size );
      if ( ( size = min( size, nb_alg ) ) == 0 ) break;

      PLA_Obj_split_4( A_TL, -size, -size, &A_TL,     &A_12,
			                   PLA_DUMMY, &A_22 );

      PLA_Obj_vert_split_2( B_L, -size, &B_L, &B_2 );
      
      PLA_Mscalar_create_conf_to( A_22, PLA_ALL_ROWS, PLA_ALL_COLS,
				  &A_22_dup );

      PLA_Obj_vert_split_2( A_12t_dup, -size, &A_12t_dup, PLA_DUMMY );

      if ( size != nb_alg ){
	PLA_Obj_vert_split_2( C_2_dup,   size, &C_2_dup, PLA_DUMMY );
	PLA_Obj_horz_split_2( A_12t_dup, size, &A_12t_dup, 
			                          PLA_DUMMY );
      }

      PLA_Copy( B_2, C_2_dup );       
      PLA_Copy( A_22, A_22_dup );

      PLA_Local_trsm( PLA_SIDE_RIGHT, PLA_UPPER_TRIANGULAR,
		       PLA_TRANSPOSE, diag,
		       one, A_22_dup, C_2_dup );
      
      PLA_Copy( A_12, A_12t_dup );

      PLA_Copy( C_2_dup, B_2 );
      
      PLA_Local_gemm( PLA_NO_TRANS, PLA_NO_TRANS,
		       minus_one, C_2_dup, A_12t_dup, one, B_L );
    }
  }

  PLA_Obj_free( &A_22 );
  PLA_Obj_free( &A_22_dup );
  PLA_Obj_free( &A_12 );
  PLA_Obj_free( &A_TL );
  PLA_Obj_free( &A_12t_dup );
  PLA_Obj_free( &B_L );
  PLA_Obj_free( &B_2 );
  PLA_Obj_free( &C_2_dup );
  PLA_Obj_free( &one );
  PLA_Obj_free( &minus_one );
  
  return PLA_SUCCESS;
}


