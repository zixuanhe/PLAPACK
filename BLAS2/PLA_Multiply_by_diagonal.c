/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Multiply_by_diagonal( int side, PLA_Obj diag, PLA_Obj B )
{
  int
    objtype;

  PLA_Obj
    diag_mv = NULL,  diag_dup = NULL;

  PLA_Obj_objtype( diag, &objtype );
  if ( objtype == PLA_MATRIX ){
    /* Copy the diagonal of the matrix to a vector */
    PLA_Mvector_create_conf_to( diag, 1, &diag_mv );

    PLA_Copy_diag_to_mv( diag, diag_mv );
  }
  else
    PLA_Obj_view_all( diag, &diag_mv );

  PLA_Obj_objtype( B, &objtype );

  if ( objtype == PLA_MVECTOR ){
    if ( side != PLA_SIDE_LEFT )
      PLA_Abort("PLA_Multiply_by_diagonal: can only multiply by diagonal from right for multivector", __LINE__, __FILE__ );

    PLA_Local_multiply_by_diagonal( side, diag_mv, B );
  }
  else if ( objtype == PLA_MATRIX ){
    if ( side == PLA_SIDE_LEFT )
      PLA_Pmvector_create_conf_to( B, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, 1, &diag_dup );
    else
      PLA_Pmvector_create_conf_to( B, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, 1, &diag_dup );

    PLA_Copy( diag_mv, diag_dup );

    PLA_Local_multiply_by_diagonal( side, diag_dup, B );
  }
  
  return PLA_SUCCESS;
}

int PLA_Local_multiply_by_diagonal( int side, PLA_Obj diag_mv, PLA_Obj B )
{
  MPI_Datatype
    datatype;

  int
    length_B, width_B, ldim_B,
    length_diag, stride_diag,
    i, j, i_one = 1;

  PLA_Obj_datatype    ( B, &datatype );
  PLA_Obj_local_length( B, &length_B );
  PLA_Obj_local_width ( B, &width_B );
  PLA_Obj_local_ldim  ( B, &ldim_B );
  
  PLA_Obj_local_length( diag_mv, &length_diag );
  if ( length_diag == 1 ) {
    PLA_Obj_local_width( diag_mv, &length_diag );
    PLA_Obj_local_ldim ( diag_mv, &stride_diag );
  }
  else{
    stride_diag = 1;
  }
      
  if ( MPI_DOUBLE == datatype ){
    double 
      *buff_B, *buff_diag;

    PLA_Obj_local_buffer( B,        ( void ** ) &buff_B );
    PLA_Obj_local_buffer ( diag_mv, ( void ** ) &buff_diag );
    if ( side == PLA_SIDE_LEFT ){
      for ( i=0; i<length_B; i++ ){
	PLA_dscal( &width_B, &buff_diag[ i*stride_diag ], &buff_B[ i ], &ldim_B );
      }
    }
    else{ /* side == PLA_SIDE_RIGHT */
      for ( j=0; j<width_B; j++ )
	PLA_dscal( &length_B, &buff_diag[ j*stride_diag ], &buff_B[ j*ldim_B ], &i_one );
    }
  }
  else if ( MPI_FLOAT == datatype ){
    float
      *buff_B, *buff_diag;

    PLA_Obj_local_buffer( B,        ( void ** ) &buff_B );
    PLA_Obj_local_buffer ( diag_mv, ( void ** ) &buff_diag );
    if ( side == PLA_SIDE_LEFT ){
      for ( i=0; i<length_B; i++ ){
	PLA_sscal( &width_B, &buff_diag[ i*stride_diag ], &buff_B[ i ], &ldim_B );
      }
    }
    else{ /* side == PLA_SIDE_RIGHT */
      for ( j=0; j<width_B; j++ )
	PLA_sscal( &length_B, &buff_diag[ j*stride_diag ], &buff_B[ j*ldim_B ], &i_one );
    }
  }
  else if ( MPI_DOUBLE_COMPLEX == datatype ){
    PLA_DOUBLE_COMPLEX
      *buff_B, *buff_diag;

    PLA_Obj_local_buffer( B,        ( void ** ) &buff_B );
    PLA_Obj_local_buffer ( diag_mv, ( void ** ) &buff_diag );
    if ( side == PLA_SIDE_LEFT ){
      for ( i=0; i<length_B; i++ ){
	PLA_zscal( &width_B, &buff_diag[ i*stride_diag ], &buff_B[ i ], &ldim_B );
      }
    }
    else{ /* side == PLA_SIDE_RIGHT */
      for ( j=0; j<width_B; j++ )
	PLA_zscal( &length_B, &buff_diag[ j*stride_diag ], &buff_B[ j*ldim_B ], &i_one );
    }
  }
  else if ( MPI_COMPLEX == datatype ){
    PLA_COMPLEX
      *buff_B, *buff_diag;

    PLA_Obj_local_buffer( B,        ( void ** ) &buff_B );
    PLA_Obj_local_buffer ( diag_mv, ( void ** ) &buff_diag );
    if ( side == PLA_SIDE_LEFT ){
      for ( i=0; i<length_B; i++ ){
	PLA_cscal( &width_B, &buff_diag[ i*stride_diag ], &buff_B[ i ], &ldim_B );
      }
    }
    else{ /* side == PLA_SIDE_RIGHT */
      for ( j=0; j<width_B; j++ )
	PLA_cscal( &length_B, &buff_diag[ j*stride_diag ], &buff_B[ j*ldim_B ], &i_one );
    }
  }

  return PLA_SUCCESS;
}
