/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Matrix_infinity_norm( PLA_Obj A, PLA_Obj alpha )
{
  int 
    size, proj_onto_A,
    value = PLA_SUCCESS;

  PLA_Obj
    A_B = NULL,      a_1 = NULL,
    row_asums = NULL, local_row_asums = NULL, local_row_asums_B = NULL,
    alpha_1 = NULL, local_max = NULL;

  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Matrix_infinity_norm_enter( A, alpha );

  /* Create a vector to hold the one norm of the columns of the matrix */

  PLA_Obj_get_orientation( A, &proj_onto_A );
  PLA_Obj_set_orientation( A, PLA_PROJ_ONTO_COL );
  PLA_Mvector_create_conf_to( A, 1, &row_asums );
  PLA_Obj_set_orientation( A, proj_onto_A );

  /* Create a vector projected onto rows to hold the one norm of the
     rows of the local matrices */

  PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1,
				&local_row_asums );

  PLA_Obj_view_all( A, &A_B );
  PLA_Obj_view_all( local_row_asums, &local_row_asums_B );

  while ( TRUE ){
    PLA_Obj_global_length( A_B, &size );
    if ( size == 0 ) break;

    PLA_Obj_horz_split_2( A_B, 1,               &a_1, 
			                        &A_B );
    PLA_Obj_horz_split_2( local_row_asums_B, 1, &alpha_1, 
                                                &local_row_asums_B );
    
    PLA_Local_asum( a_1, alpha_1 );
  }
  
  /* compute global 1-norms of the rows */
  PLA_Reduce( local_row_asums, MPI_SUM, row_asums );
  
  /* Create a multiscalar in which to hold the local asum */

  PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS, 
			       &local_max );

  PLA_Local_absolute_max( row_asums, local_max );

  PLA_Reduce( local_max, MPI_MAX, alpha );

  PLA_Obj_free ( &A_B );
  PLA_Obj_free ( &a_1 );
  PLA_Obj_free ( &row_asums );
  PLA_Obj_free ( &local_row_asums );
  PLA_Obj_free ( &local_row_asums_B );
  PLA_Obj_free ( &local_max );
  PLA_Obj_free ( &alpha_1 );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Matrix_infinity_norm_exit( A, alpha );

}
