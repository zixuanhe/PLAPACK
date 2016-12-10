/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_bidiag_svd( PLA_Obj bidiag, PLA_Obj U_mv, PLA_Obj V_mv )
{
  int 
    info, length_D,
    local_length_U, local_length_V,
    i_zero = 0, i_one = 1, ldU, ldV;
  void 
    *buffer_diag, *buffer_superdiag, *buffer_U, *buffer_V, *work;
  MPI_Datatype
    datatype;
  PLA_Obj
    diag = NULL, superdiag = NULL;

  PLA_Obj_vert_split_2( bidiag, 1, &diag, &superdiag );

  PLA_Obj_datatype( diag, &datatype );
  PLA_Obj_global_length( diag, &length_D );
  PLA_Obj_local_buffer( diag,      &buffer_diag );
  PLA_Obj_local_buffer( superdiag, &buffer_superdiag );

  if ( U_mv == NULL )
    local_length_U = 0;
  else{
    PLA_Obj_local_length( U_mv, &local_length_U );
    PLA_Obj_local_buffer( U_mv, &buffer_U );
    PLA_Obj_local_ldim  ( U_mv, &ldU );
  }

  if ( V_mv == NULL )
    local_length_V = 0;
  else{
    PLA_Obj_local_length( V_mv, &local_length_V );
    PLA_Obj_local_buffer( V_mv, &buffer_V );
    PLA_Obj_local_ldim  ( V_mv, &ldV );
  }

  work = ( double * ) PLA_malloc( 4 * length_D * sizeof( double ) );
  PLA_dbdsqr_x( "Upper", &length_D, &local_length_V, &local_length_U, &i_zero, 
	      buffer_diag, buffer_superdiag, buffer_V, &ldV, buffer_U, &ldU,
	      NULL, &i_one, work, &info );

  if ( info < 0 )
    PLA_Abort( "dbdsqr_x returned with parameter error", __LINE__, __FILE__ );
  else if ( info > 0 )
    PLA_Abort( "dbdsqr_x returned convergence problem", __LINE__, __FILE__ );

  PLA_free( work );
  PLA_Obj_free( &diag );
  PLA_Obj_free( &superdiag );

  return PLA_SUCCESS;
}
