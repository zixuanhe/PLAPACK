/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Gemm( int transa, int transb, 
	       PLA_Obj alpha, PLA_Obj A, PLA_Obj B, 
               PLA_Obj beta,  PLA_Obj C )

/*
  Purpose : Parallel matrix multiplication

  IN     transa      integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     transb      integer, PLA_NO_TRANSPOSE, PLA_TRANSPOSE, 
                              or PLA_CONJUGATE TRANSPOSE
  IN     alpha       multiscalar, scale factor for A and B
  IN     A           matrix
  IN     B           matrix
  IN     beta        multiscalar, scale factor for C
  IN/OUT C           matrix, overwritten with 
                     C <- alpha * A  * B  + beta * C  or 
                     C <- alpha * A  * B' + beta * C  or 
                     C <- alpha * A' * B  + beta * C  or 
                     C <- alpha * A  * B' + beta * C

  NOTE:  For details on how to implement matrix-matrix multiplication,
         see
	 
	 R. van de Geijn, Using PLAPACK, The MIT Press, 1997.

	 R. van de Geijn and J. Watts, 
	 "SUMMA: Scalable Universal Matrix Multiplication Algorith,"
	 Concurrency: Practice and Experience, Vol 9 (4), pp. 255-274
	 (April 1997).
	 
	 G. Morrow, J. Gunnels, C. Lin, and R. van de Geijn,
	 "A Flexible Class of Parallel Matrix Multiplication Algorithms,"
	 Proceedings of IPPS98, pp. 110-116, 1998.
*/
{
  int 
    value = PLA_SUCCESS,
    owner_row, owner_col,
    length_A, width_A, 
    length_B, width_B, 
    length_C, width_C, 
    nb_alg;
  
  PLA_Obj
    alpha_cpy = NULL, beta_cpy = NULL;

  PLA_Template
    templ = NULL;

  /* Perform parameter and error checking */
  if ( PLA_ERROR_CHECKING )    
    value = PLA_Gemm_enter( transa, transb, alpha, A, B, beta, C );

  if ( !value ){
    /* If necessary, duplicate alpha and beta to all nodes */
    PLA_Obj_owner_row( alpha, &owner_row );
    PLA_Obj_owner_col( alpha, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &alpha_cpy );
      PLA_Copy( alpha, alpha_cpy );
    }

    PLA_Obj_owner_row( beta, &owner_row );
    PLA_Obj_owner_col( beta, &owner_col );
    if ( owner_row != PLA_ALL_ROWS || owner_col != PLA_ALL_COLS ){
      PLA_Mscalar_create_conf_to( beta, PLA_ALL_ROWS, PLA_ALL_COLS,
				   &beta_cpy );
      PLA_Copy( beta, beta_cpy );
    }

    /* Get template to be used to extract algorithmic block size, later */
    PLA_Obj_template( A, &templ );

    /* Extract the dimensions of the different matrices */
    PLA_Obj_global_length( A, &length_A );
    PLA_Obj_global_width(  A, &width_A );
    PLA_Obj_global_length( B, &length_B );
    PLA_Obj_global_width(  B, &width_B );
    PLA_Obj_global_length( C, &length_C );
    PLA_Obj_global_width(  C, &width_C );

    if ( length_A * width_A > length_B * width_B && 
	 length_A * width_A > length_C * width_C ){
      /* Matrix A has most data.  Leave it in place, communicating
         the other two */

      /* Get algorithmic block size to be used */
      PLA_Environ_nb_alg( PLA_OP_MAT_PAN, templ, &nb_alg );

      PLA_Gemm_A( nb_alg, transa, transb, 
		      ( alpha_cpy == NULL ? alpha: alpha_cpy ), A, B, 
		      ( beta_cpy  == NULL ? beta : beta_cpy  ), C ); 
    }
    else if ( length_B * width_B > length_A * width_A && 
	      length_B * width_B > length_C * width_C ){
      /* Matrix B has most data.  Leave it in place, communicating
         the other two */

      /* Get algorithmic block size to be used */
      PLA_Environ_nb_alg( PLA_OP_PAN_MAT, templ, &nb_alg );

      PLA_Gemm_B( nb_alg, transa, transb, 
		      ( alpha_cpy == NULL ? alpha: alpha_cpy ), A, B, 
		      ( beta_cpy  == NULL ? beta : beta_cpy  ), C ); 
    }
    else {
      /* Matrix C has most data.  Leave it in place, communicating
         the other two */

      /* Get algorithmic block size to be used */
      PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

      PLA_Gemm_C( nb_alg, transa, transb, 
		      ( alpha_cpy == NULL ? alpha: alpha_cpy ), A, B, 
		      ( beta_cpy  == NULL ? beta : beta_cpy  ), C ); 
    }

    /* free temporary objects */
    PLA_Obj_free( &alpha_cpy );
    PLA_Obj_free( &beta_cpy );
  }

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Gemm_exit( transa, transb, alpha, A, B, beta, C );

  return value;
}


