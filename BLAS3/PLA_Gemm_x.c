/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Gemm_x( int version, int transa, int transb, 
	         PLA_Obj alpha, PLA_Obj A, PLA_Obj B, 
                 PLA_Obj beta,  PLA_Obj C )
/*
  Purpose : Parallel matrix multiplication 

  IN     version     version to be used.
  IN     transa      integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     transb      integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     alpha       multiscalar, scale factor for A and B
  IN     A           matrix, factor in product
  IN     B           matrix, factor in product
  IN     beta        multiscalar, scale factor for C
  IN/OUT C           matrix, overwritten with 
                     C <- alpha * A  * B  + beta * C  or 
                     C <- alpha * A  * B' + beta * C  or 
                     C <- alpha * A' * B  + beta * C  or 
                     C <- alpha * A  * B' + beta * C

*/
{
  int 
    value = PLA_SUCCESS,
    owner_row, owner_col,
    nb_alg;
  
  PLA_Obj
    alpha_cpy = NULL, beta_cpy = NULL;

  PLA_Template
    templ = NULL;


  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Gemm_x_enter( version, transa, transb, alpha, A, B, beta, C );


  if ( value == PLA_SUCCESS ){
    PLA_Obj_template( A, &templ );

    /* If necessary, duplicate alpha and beta to all nodes */
    PLA_Obj_owner_row( alpha, &owner_row );
    PLA_Obj_owner_col( alpha, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &alpha_cpy );
      PLA_Copy( beta, beta_cpy );
    }

    PLA_Obj_owner_row( beta, &owner_row );
    PLA_Obj_owner_col( beta, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( beta, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &beta_cpy );
      PLA_Copy( beta, beta_cpy );
    }

    if ( version == PLA_VERSION_C ){
      /* Matrix C has most data.  Leave it in place, communicating
         the other two */

      PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );
      PLA_Gemm_C( nb_alg, transa, transb, 
		      ( alpha_cpy == NULL ? alpha: alpha_cpy ), A, B, 
		      ( beta_cpy  == NULL ? beta : beta_cpy  ), C );
    }
    else if ( version == PLA_VERSION_A ){
      /* Matrix A has most data.  Leave it in place, communicating
         the other two */

      PLA_Environ_nb_alg( PLA_OP_MAT_PAN, templ, &nb_alg );
      PLA_Gemm_A( nb_alg, transa, transb, 
		      ( alpha_cpy == NULL ? alpha: alpha_cpy ), A, B, 
		      ( beta_cpy  == NULL ? beta : beta_cpy  ), C );
    }
    else if ( version == PLA_VERSION_B ){
      /* Matrix B has most data.  Leave it in place, communicating
         the other two */

      PLA_Environ_nb_alg( PLA_OP_PAN_MAT, templ, &nb_alg );
      PLA_Gemm_B( nb_alg, transa, transb, 
		      ( alpha_cpy == NULL ? alpha: alpha_cpy ), A, B, 
		      ( beta_cpy  == NULL ? beta : beta_cpy  ), C );
    }
    else
      PLA_Abort( "Type not supported by PLA_Gemm_x", __LINE__, __FILE__ );

    /* free temporary objects */
    PLA_Obj_free( &alpha_cpy );
    PLA_Obj_free( &beta_cpy );
  }

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Gemm_x_exit( version, transa, transb, alpha, A, B, beta, C );

  return value;
}

