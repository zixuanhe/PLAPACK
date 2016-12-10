/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Apply_W_Y_transform (int side, int trans,
			     PLA_Obj W, PLA_Obj Y, PLA_Obj A)
     /***************************************************
       OPERATION             SIDE               TRANS
       ----------------------------------------------
       A <-- (1 + W Y') A    PLA_SIDE_LEFT      PLA_NO_TRANS
       A <-- A (1 + W Y')    PLA_SIDE_RIGHT     PLA_NO_TRANS
       A <-- (1 + Y W') A    PLA_SIDE_LEFT      PLA_TRANS
       A <-- A (1 + Y W')    PLA_SIDE_RIGHT     PLA_TRANS
     ***************************************************/
{
  int
    global_width_A = 0,
    global_length_A = 0;

  PLA_Obj_global_length( A, &global_length_A);
  PLA_Obj_global_width( A, &global_width_A);

  if ( global_length_A*global_width_A == 0)
    return PLA_SUCCESS;

  if ( PLA_TRANSPOSE == trans )
    {
      if ( PLA_SIDE_LEFT == side )
	{
	  return pla_Apply_W_Y_transform_L_N( Y, W, A);
	}
      else if ( PLA_SIDE_RIGHT == side )
	{
	  return pla_Apply_W_Y_transform_R_N( Y, W, A);
	}
      else
	{
	  PLA_Abort("PLA_Apply_W_Y_transform: Invalid side specification", __LINE__, __FILE__ );
	}
    }
  else if ( PLA_NO_TRANSPOSE == trans )
    {
      if ( PLA_SIDE_LEFT == side )
	{
	  return pla_Apply_W_Y_transform_L_N ( W, Y, A);
	}
      else if ( PLA_SIDE_RIGHT == side )
	{
	  return pla_Apply_W_Y_transform_R_N( W, Y, A);
	}
      else
	{
	  PLA_Abort( "PLA_Apply_W_Y_transform: Invalid side specification", __LINE__, __FILE__ );
	}
    }
  else
    {
      PLA_Abort( "PLA_Apply_W_Y_transform: Invalid transpose specification", __LINE__, __FILE__ );
    }

}


int pla_Apply_W_Y_transform_L_N ( PLA_Obj W, PLA_Obj Y, PLA_Obj A)
{
  int 
    width;

  PLA_Obj
    W_dpmv = NULL,
    Y_dpmv = NULL,
    U_dpmv = NULL,
    U_dpmv_temp = NULL,
    one = NULL,
    zero = NULL;

  PLA_Create_constants_conf_to( W, NULL, &zero, &one );

  PLA_Obj_global_width( W, &width );

  PLA_Pmvector_create_conf_to (A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
			       width, &W_dpmv);
			       
  PLA_Pmvector_create_conf_to (A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
			       width, &Y_dpmv);

  PLA_Pmvector_create_conf_to (A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
			       width, &U_dpmv);

  PLA_Pmvector_create_conf_to (A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
			       width, &U_dpmv_temp);

  PLA_Copy ( W, W_dpmv);
  PLA_Copy ( Y, Y_dpmv);

  PLA_Local_gemm ( PLA_TRANS, PLA_NO_TRANS, 
		   one, Y_dpmv, A, zero, U_dpmv_temp);

  PLA_Reduce ( U_dpmv_temp, MPI_SUM, U_dpmv);

  PLA_Local_gemm ( PLA_NO_TRANS, PLA_NO_TRANS, 
		   one, W_dpmv, U_dpmv, one, A);

  PLA_Obj_free ( &W_dpmv);
  PLA_Obj_free ( &Y_dpmv);
  PLA_Obj_free ( &U_dpmv);
  PLA_Obj_free ( &U_dpmv_temp);
  PLA_Obj_free ( &one);
  PLA_Obj_free ( &zero);

  return PLA_SUCCESS;
}



int pla_Apply_W_Y_transform_R_N ( PLA_Obj W, PLA_Obj Y, PLA_Obj A)
{
  int
    width;

  PLA_Obj
    W_dpmv = NULL,
    Y_dpmv = NULL,
    U_dpmv = NULL,
    U_dpmv_temp = NULL,
    one = NULL,
    zero = NULL;

  PLA_Create_constants_conf_to( W, NULL, &zero, &one );

  PLA_Obj_global_width( W, &width );

  PLA_Pmvector_create_conf_to (A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
			       width, &W_dpmv);
			       
  PLA_Pmvector_create_conf_to (A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
			       width, &Y_dpmv);

  PLA_Pmvector_create_conf_to (A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
			       width, &U_dpmv);

  PLA_Pmvector_create_conf_to (A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
			       width, &U_dpmv_temp);

  PLA_Copy ( W, W_dpmv);
  PLA_Copy ( Y, Y_dpmv);

  PLA_Local_gemm ( PLA_NO_TRANS, PLA_TRANS, 
		   one, A, W_dpmv, zero, U_dpmv_temp);

  PLA_Reduce ( U_dpmv_temp, MPI_SUM, U_dpmv);

  PLA_Local_gemm ( PLA_NO_TRANS, PLA_NO_TRANS, 
		   one, U_dpmv, Y_dpmv, one, A);

  PLA_Obj_free ( &W_dpmv);
  PLA_Obj_free ( &Y_dpmv);
  PLA_Obj_free ( &U_dpmv);
  PLA_Obj_free ( &U_dpmv_temp);
  PLA_Obj_free ( &one);
  PLA_Obj_free ( &zero);

  return PLA_SUCCESS;
}

