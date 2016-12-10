/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Matrix_one_norm( PLA_Obj A, PLA_Obj alpha )
{
  int 
    size, proj_onto_A,
    value = PLA_SUCCESS;

  PLA_Obj
    A_R = NULL,      a_1 = NULL,
    col_asums = NULL, local_col_asums = NULL, local_col_asums_R = NULL,
    alpha_1 = NULL, local_max = NULL;

  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Matrix_one_norm_enter( A, alpha );

  /* Create a vector to hold the one norm of the columns of the matrix */

  PLA_Obj_get_orientation( A, &proj_onto_A );
  PLA_Obj_set_orientation( A, PLA_PROJ_ONTO_ROW );
  PLA_Mvector_create_conf_to( A, 1, &col_asums );
  PLA_Obj_set_orientation( A, proj_onto_A );

  /* Create a vector projected onto rows to hold the one norm of the
     columns of the local matrices */

  PLA_Pmvector_create_conf_to( A, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1,
				&local_col_asums );

  PLA_Obj_view_all( A, &A_R );
  PLA_Obj_view_all( local_col_asums, &local_col_asums_R );

  while ( TRUE ){
    PLA_Obj_global_width( A_R, &size );
    if ( size == 0 ) break;

    PLA_Obj_vert_split_2( A_R, 1,               &a_1,     &A_R );
    PLA_Obj_vert_split_2( local_col_asums_R, 1,  &alpha_1, &local_col_asums_R );
    
    PLA_Local_asum( a_1, alpha_1 );
  }
  
  /* compute global 1-norms of the columns */
  PLA_Reduce( local_col_asums, MPI_SUM, col_asums );
  
  /* Create a multiscalar in which to hold the local asum */

  PLA_Mscalar_create_conf_to( alpha, PLA_ALL_ROWS, PLA_ALL_COLS, 
			       &local_max );

  PLA_Local_absolute_max( col_asums, local_max );

  PLA_Reduce( local_max, MPI_MAX, alpha );

  PLA_Obj_free ( &A_R );
  PLA_Obj_free ( &a_1 );
  PLA_Obj_free ( &col_asums );
  PLA_Obj_free ( &local_col_asums );
  PLA_Obj_free ( &local_col_asums_R );
  PLA_Obj_free ( &local_max );
  PLA_Obj_free ( &alpha_1 );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Matrix_one_norm_exit( A, alpha );

}
