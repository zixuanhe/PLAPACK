/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Nrm2( PLA_Obj x, PLA_Obj alpha )
{
  int 
    nprocs, owner_row, owner_col,
    value = PLA_SUCCESS;

  void 
    *buff_all, *buff_local;

  PLA_Obj
    alpha_all = NULL, alpha_local = NULL;

  MPI_Datatype
    datatype;

  MPI_Comm
    comm = MPI_COMM_NULL;

  PLA_Template
    templ;

  if ( PLA_ERROR_CHECKING )    
    value = PLA_Nrm2_enter( x, alpha );

  PLA_Obj_template( alpha, &templ );
  PLA_Temp_comm_all_size( templ, &nprocs );
  PLA_Temp_comm_all( templ, &comm );

  PLA_Obj_datatype( alpha, &datatype );
  PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS, 
		       nprocs, 1, templ, &alpha_all );

  PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS, 
		       1, 1, templ, &alpha_local );

  PLA_Local_nrm2( x, alpha_local );
  PLA_Obj_local_buffer( alpha_local, &buff_local );
  PLA_Obj_local_buffer( alpha_all, &buff_all );
    
  MPI_Allgather( buff_local, 1, datatype, buff_all, 1, datatype,
		 comm );

  PLA_Obj_owner_row( alpha, &owner_row );
  PLA_Obj_owner_col( alpha, &owner_col );

  if ( owner_row == PLA_ALL_ROWS && owner_col == PLA_ALL_COLS )
    PLA_Local_nrm2( alpha_all, alpha );
  else{
    PLA_Local_nrm2( alpha_all, alpha_local );
    PLA_Copy( alpha_local, alpha );
  }

  PLA_Obj_free( &alpha_all );
  PLA_Obj_free( &alpha_local );

  if ( PLA_ERROR_CHECKING )   
    value = PLA_Nrm2_exit( x, alpha );

  return value;
}


