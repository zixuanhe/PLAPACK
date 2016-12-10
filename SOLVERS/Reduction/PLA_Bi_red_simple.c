/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Bi_red( int uplo, PLA_Obj A, PLA_Obj sL, PLA_Obj sR, PLA_Obj U, PLA_Obj V )
/*
  PLA_Tri_red

  Purpose: Reduce general matrix A to bidiagonal form using 
  Householder similarity transformations.

  input:
  A                       MATRIX to be reduced

  output:
  A                       Reduced matrix A.  Householder vectors used
                          to reduce A are stored below subdiagonal
                          of A and above first superdiagonal of A.
  sL                      Scaling factors for the Householder transforms
                          computed to reduce A (applied to left of A).
  sR                      Scaling factors for the Householder transforms
                          computed to reduce A (applied to right of A).
  U                       if U != NULL, U equals the accumulation of 
                          Householder transforms that were applied from the
                          left.
  V                       if V != NULL, V equals the accumulation of 
                          Householder transforms that were applied from the
                          right.
*/
{
  PLA_Obj
    uL = NULL,  uL_B = NULL, uL_1 = NULL, 
    uR = NULL,  uR_B = NULL, uR_1 = NULL, 
    sR_B = NULL, betaR_1 = NULL, betaR_B = NULL,  
    sL_B = NULL, betaL_1 = NULL, betaL_B = NULL,  
    beta_1_dup = NULL,
    a_L = NULL, A_R = NULL, A_BR = NULL, a_12 = NULL; 
  int
    length, width, size, value = PLA_SUCCESS;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Bi_red_enter( uplo, A, sL, sR, U, V );

  if ( uplo == PLA_UPPER_TRIANGULAR ){
    PLA_Mvector_create_conf_to( A, 1, &uL );
    PLA_Obj_set_orientation( A, PLA_PROJ_ONTO_ROW );
    PLA_Mvector_create_conf_to( A, 1, &uR );
    PLA_Obj_set_orientation( A, PLA_PROJ_ONTO_COL );
  
    /* Create duplicated multiscalars in which to hold the scaling factor
       for the Householder transforms being computed */
    PLA_Obj_horz_split_2( sL, 1, &betaL_1,
			         PLA_DUMMY );
    PLA_Mscalar_create_conf_to( betaL_1, PLA_ALL_ROWS, PLA_ALL_COLS, 
			         &beta_1_dup );

    /* Track the active parts of A, sL, sR, and uwR and wLv */
    PLA_Obj_view_all( A, &A_BR );    
    PLA_Obj_view_all( sL, &betaL_B );    

    PLA_Obj_horz_split_2( sR, 0, PLA_DUMMY,
			         &betaR_B );

    PLA_Obj_view_all( uL, &uL_B );
    PLA_Obj_view_all( uR, &uR_B );

    while ( TRUE ){  
      PLA_Obj_global_width( A_BR, &size );

      if ( 0 == size ) break;

      /* Partition A_BR = ( a_L      A_R    ) where a_L is a column from
	 which to compute the next u */
      PLA_Obj_vert_split_2( A_BR, 1, &a_L, &A_R );
      /* Split of the current element of vector sL */
      PLA_Obj_horz_split_2( betaL_B, 1, &betaL_1,
                                        &betaL_B );

      /* Split of the current element of vector sR */
      PLA_Obj_horz_split_2( betaR_B, 1, &betaR_1,
                                        &betaR_B );

      /* Redistributed a_L as a vector and compute Householder transform */
      PLA_Copy( a_L, uL_B );  
      PLA_Compute_House_v( uL_B, beta_1_dup );

      /* Place data back in A and sL */
      PLA_Copy( uL_B, a_L );
      PLA_Local_copy( beta_1_dup, betaL_1 );

      /* Set first element of uL_B to one */
      PLA_Obj_horz_split_2( uL_B, 1, &uL_1,
			             PLA_DUMMY );
      PLA_Obj_set_to_one( uL_1 );

      /* Apply Householder transform to A_R  from left */
      PLA_Apply_House_transform( PLA_SIDE_LEFT, A_R, uL_B, beta_1_dup );

      PLA_Obj_horz_split_2( uL_B, 1, PLA_DUMMY,
			             &uL_B );

      if ( 1 == size ) break;
    
      /* Partition A_BR = /  *    a_12^T \  
                          \  *    A_22   / splitting off first column and row */

      PLA_Obj_split_4( A_BR, 1, 1, PLA_DUMMY, &a_12,
	   	                   PLA_DUMMY, &A_BR );

      /* Update u_B so we can copy a_12 into it */

      PLA_Obj_horz_split_2( uR_B, 1, PLA_DUMMY,
			             &uR_B );

      /* Copy a_12 into uR_B */
      PLA_Obj_set_orientation( a_12, PLA_PROJ_ONTO_ROW );
      PLA_Copy( a_12, uR_B );

      /* Compute Householder transform */
      PLA_Compute_House_v( uR_B, beta_1_dup );

      /* Place data back in A and sL */
      PLA_Copy( uR_B, a_12 ); 
      PLA_Local_copy( beta_1_dup, betaR_1 ); 

      /* Set first element of uR_B to one */
      PLA_Obj_horz_split_2( uR_B, 1, &uR_1,
			             PLA_DUMMY );

      PLA_Obj_set_to_one( uR_1 );

      /* Apply Householder transform to A_BR  from right */
      PLA_Apply_House_transform( PLA_SIDE_RIGHT, A_BR, uR_B, beta_1_dup );
    }
  }
  else 
    PLA_Abort( "PLA_LOWER_TRIANGULAR not yet supported",
	       __LINE__, __FILE__ );

  /* Compute U and V */
  if ( U != NULL )
    PLA_Form_Q( PLA_NO_TRANSPOSE, A, sL, U );

  if ( V != NULL ){
    PLA_Obj
      A12 = NULL, V22 = NULL;

    PLA_Obj_vert_split_2( A, 1, PLA_DUMMY, &A12 );

    PLA_Obj_set_to_identity( V );
    PLA_Obj_split_4( V, 1, 1, PLA_DUMMY, PLA_DUMMY,
		              PLA_DUMMY, &V22 );

    PLA_Form_Q_x( PLA_NO_TRANSPOSE, PLA_UPPER_TRIANGULAR,
		  A12, sR, V22 );

    PLA_Obj_free( &A12 );
    PLA_Obj_free( &V22 );
  }

  /* Free the temporary objects */

  PLA_Obj_free( &uL );
  PLA_Obj_free( &uL_B );
  PLA_Obj_free( &uL_1 );
  PLA_Obj_free( &uR );
  PLA_Obj_free( &uR_B );
  PLA_Obj_free( &uR_1 );
  PLA_Obj_free( &sR_B );
  PLA_Obj_free( &betaR_1 );
  PLA_Obj_free( &betaR_B );
  PLA_Obj_free( &sL_B );
  PLA_Obj_free( &betaL_1 );
  PLA_Obj_free( &betaL_B );
  PLA_Obj_free( &beta_1_dup );
  PLA_Obj_free( &a_L );
  PLA_Obj_free( &A_R );
  PLA_Obj_free( &A_BR );
  PLA_Obj_free( &a_12 );


  if ( PLA_ERROR_CHECKING )   
    value = PLA_Bi_red_exit( uplo, A, sL, sR, U, V );

  return PLA_SUCCESS;
}
