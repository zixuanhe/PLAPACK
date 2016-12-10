/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Local_sym_tridiag_eig( PLA_Obj tridiag, 
				PLA_Obj Z )
{
  int 
    info, length, i_one = 1, local_length_Z, ldz;
  void 
    *buffer_diag, *buffer_subdiag, *buffer_Z, *work;
  MPI_Datatype
    datatype;
  PLA_Obj
    diag = NULL, subdiag = NULL;

  PLA_Obj_vert_split_2( tridiag, 1, &diag, &subdiag );

  PLA_Obj_datatype( diag, &datatype );
  PLA_Obj_global_length( diag, &length );
  PLA_Obj_local_buffer( diag,    &buffer_diag );
  PLA_Obj_local_buffer( subdiag, &buffer_subdiag );

  if ( Z == NULL )
    PLA_dsteqr( "N", &length, 
		   (double *) buffer_diag, (double *) buffer_subdiag,
		   NULL, &i_one, work, &info );
  else{
    PLA_Obj_local_length( Z, &local_length_Z );
    PLA_Obj_local_buffer( Z, &buffer_Z );
    PLA_Obj_local_ldim( Z, &ldz );
    work = ( double * ) PLA_malloc( 2 * length * sizeof( double ) );
    PLA_dsteqr_x( "V", &local_length_Z, &length, 
		   (double *) buffer_diag, (double *) buffer_subdiag,
		   (double *) buffer_Z, &ldz, work, &info );
    PLA_free( work );
  }

  PLA_Obj_free( &diag );
  PLA_Obj_free( &subdiag );

  return PLA_SUCCESS;
}
