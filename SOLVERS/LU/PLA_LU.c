/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_LU( PLA_Obj A, PLA_Obj pivots )
{
  int 
    value = PLA_SUCCESS,
    nb_alg;

  PLA_Template
    templ;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_LU_enter( A, pivots );

  PLA_Obj_template( A, &templ );
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  value = PLA_LU_right( nb_alg, A, pivots );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_LU_exit( A, pivots );


  return value;
}

int PLA_LU_right( int nb_alg, PLA_Obj A, PLA_Obj pivots )
/*
   Blocked right-looking LU factorization;

   Input: 

   nb_alg        --  Algorithmic blocking size to be used
  
   A             --  m x n PLA_Matrix
                     contains original matrix to be factored

   pivots        --  PLA_Mvector of width 1 and datatype MPI_INT
                     vector of length min( m, n ) in which to store
                     the pivot information

   Output:

   A             --  m x n PLA_Matrix
                     contains factors L and U s.t. P A_orig = L U

   pivots        --  PLA_Mvector of width 1 and datatype MPI_INT
                     integer vector that describes permutation matrix P
		     
   returns value PLA_SUCCESS if factorization completes successfully
*/
{
  int       
    value = PLA_SUCCESS,
    size, length, width;

  PLA_Obj  
    A_BR            = NULL,     A_1            = NULL,
    A_B             = NULL,
    A_11            = NULL,     A_12            = NULL,     
    A_21            = NULL,     A_21_dpmv       = NULL,
    A_1_mv          = NULL,     A_1_dpmv        = NULL,
    A_11_msc        = NULL,     A_12_dpmv       = NULL,
    A_12_mv         = NULL,
    pivots_B        = NULL,     pivots_1       = NULL,
    pivots_1_msc    = NULL,
    one             = NULL,     minus_one      = NULL;

  /* Create usual duplicated scalar constants */
  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  /* Partition A = / A_TL  A_TR \
                   \ A_BL  A_BR / where A_TL is 0 x 0 

     A_BR will track the active part of matrix A       */
  PLA_Obj_view_all( A, &A_BR );

  /* Partition A = / A_T \
                   \ A_B /   where A_T is 0 x n

     A_B will track the part of A that needs to be pivoted */

  PLA_Obj_view_all( A, &A_B );
 
  /* Partition pivots = / pivots_T \
                        \ pivots_B /   where pivots_T is 0 x 1 
    
     pivots_B will track the active part of the pivot vector */

  PLA_Obj_view_all( pivots, &pivots_B );

  /* Do until done */
  while ( TRUE ) {          
    /* Determine size of current panel and check if done */
    PLA_Obj_global_length( A_BR, &length );
    PLA_Obj_global_width ( A_BR,  &width );
    size = min( length, width );
    if ( 0 == ( size = min( size, nb_alg ) ) ) break;

    /* Partition off the current column panel to be factored */
    PLA_Obj_vert_split_2( A_BR, size,       &A_1, PLA_DUMMY );

    /* Create a multivector to hold current column panel, allowing
       all nodes to participate in the panel factorization */
    PLA_Mvector_create_conf_to( A_1, size, &A_1_mv );

    /* Copy the current column panel to the multivector */
    PLA_Copy( A_1, A_1_mv );

    /* Partition A_BR = / A_11 A_12 \
                        \ A_21 A_22 /  where A_11 is size x size */	       
    PLA_Obj_split_4( A_BR, size, size,      &A_11, &A_12,
                                             &A_21, &A_BR );

    /* Create a duplicated multiscalar to hold a copy of A_11 after
       the panel has been factored (side-effect of the panel factorization
       routine */
    PLA_Mscalar_create_conf_to( A_11, PLA_ALL_ROWS, PLA_ALL_COLS, 
                                 &A_11_msc );

    /* Partition off the subvector of pivots that will hold the pivot
       information for the column panel to be factored */
    PLA_Obj_horz_split_2( pivots_B, size,   &pivots_1, 
                                             &pivots_B );

    /* Create a duplicated multiscalar in which to store the pivot
       information for the column panel factorization */
    PLA_Mscalar_create_conf_to( pivots_1, PLA_ALL_ROWS, PLA_ALL_COLS, 
				 &pivots_1_msc );

    /* Factor the current column panel, distributed as a multivector.
       the updated A_11 matrix will be returned duplicated to all
       nodes in A_11_msc */ 	       
    value = PLA_LU_right_mv( A_1_mv, A_11_msc, pivots_1_msc  );
    if ( value != PLA_SUCCESS )  break;

    /* Copy the pivot information for the current column panel into pivots */
    PLA_Copy( pivots_1_msc, pivots_1 );

    /* Pivot the rows of the matrix consistent with the pivot information
       for the current column panel */
    PLA_Apply_pivots_to_rows( A_B, pivots_1_msc );

    /* Copy A_12 so it is duplicated to all rows, in A_12_dpmv */
    PLA_Obj_set_orientation( A_12, PLA_PROJ_ONTO_ROW );
    PLA_Pmvector_create_conf_to( A_12, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
				  size, &A_12_dpmv );

    /* Duplicate the current panel to all columns, in A_1_dpmv */
    PLA_Pmvector_create_conf_to( A_1, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
				  size, &A_1_dpmv );
    PLA_Copy( A_1_mv,   A_1_dpmv );

    /* Place the current column panel back into the matrix */
    PLA_Copy( A_1_dpmv, A_1 );  

    /* View view the part of A_1_dpmv that holds the copy of the updated
       matrix A_21 */
    PLA_Obj_horz_split_2( A_1_dpmv, size, PLA_DUMMY,
			                   &A_21_dpmv );
    
    /* Update A_12 <- inv(L_11) A_12, while distributed in A_12_dpmv */

    PLA_Mvector_create_conf_to( A_12_dpmv, size, &A_12_mv );
    PLA_Copy( A_12, A_12_mv );
    PLA_Local_trsm( PLA_SIDE_RIGHT, PLA_LOWER_TRIANGULAR, 
		     PLA_TRANSPOSE, PLA_UNIT_DIAG,
		     one, A_11_msc, A_12_mv );
    PLA_Copy( A_12_mv, A_12_dpmv ); 

    /* Place the current row panel A_12 back into the matrix */
    PLA_Copy( A_12_dpmv, A_12 );

    /* Update A_22 (viewed by A_BR) <- A_22 - A_21 * A_12 */
    PLA_Local_gemm( PLA_NO_TRANSPOSE, PLA_NO_TRANSPOSE, 
		     minus_one, A_21_dpmv, A_12_dpmv, one, A_BR ); 

    /* Update the view of the part of the matrix that still needs to
       be pivoted in subsequent steps */
    PLA_Obj_horz_split_2( A_B, size,       PLA_DUMMY,
			                    &A_B );
  }
  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &A_1 );
  PLA_Obj_free( &A_B );
  PLA_Obj_free( &A_11 );
  PLA_Obj_free( &A_12 );
  PLA_Obj_free( &A_21 );
  PLA_Obj_free( &A_21_dpmv );
  PLA_Obj_free( &A_1_mv );
  PLA_Obj_free( &A_1_dpmv );
  PLA_Obj_free( &A_12_dpmv );
  PLA_Obj_free( &A_12_mv );
  PLA_Obj_free( &A_11_msc );
  PLA_Obj_free( &pivots_B );
  PLA_Obj_free( &pivots_1 );
  PLA_Obj_free( &pivots_1_msc );
  PLA_Obj_free( &one );
  PLA_Obj_free( &minus_one );

  return value;
}

