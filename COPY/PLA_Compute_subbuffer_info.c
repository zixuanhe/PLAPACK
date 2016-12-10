/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Compute_subbuffer_info_within_row( 
       PLA_Obj Obj_mv,
       int *counts,
       int *displs )
{
  int
    align_row_mv, 
    length_mv, width_mv, 
     me, i,
     nprows, npcols, 
     myrow_first_rank; 
  PLA_Template templ = NULL;
  MPI_Datatype datatype;

  PLA_Obj_global_align_row( Obj_mv, &align_row_mv );
  PLA_Obj_global_length( Obj_mv, &length_mv );
  PLA_Obj_global_width( Obj_mv, &width_mv );

  PLA_Obj_template( Obj_mv, &templ );
  PLA_Temp_comm_all_rank( templ, &me );
  PLA_Temp_comm_row_size( templ, &npcols );
  PLA_Temp_comm_col_size( templ, &nprows );

  myrow_first_rank = me%nprows;

  displs[0] = 0;
  for ( i=0; i<npcols; i++ ) {
    PLA_Temp_vector_distant_length(templ,
                              myrow_first_rank+i*nprows, 
                              length_mv,
                              align_row_mv,
                              &counts[i]);
    counts[i] *= width_mv;
    if ( i>0 ) displs[i] = displs[i-1] + counts[i-1];
  }

  return PLA_SUCCESS;
}


int PLA_Compute_subbuffer_info_within_col( 
       PLA_Obj Obj_mv,
       int *counts,
       int *displs )
{
  int
    align_row_mv, 
    length_mv, width_mv, 
     me, i,
     nprows, npcols, 
     mycol_first_rank; 
  PLA_Template templ = NULL;
  MPI_Datatype datatype;

  PLA_Obj_global_align_row( Obj_mv, &align_row_mv );
  PLA_Obj_global_length( Obj_mv, &length_mv );
  PLA_Obj_global_width( Obj_mv, &width_mv );

  PLA_Obj_template( Obj_mv, &templ );
  PLA_Temp_comm_all_rank( templ, &me );
  PLA_Temp_comm_row_size( templ, &npcols );
  PLA_Temp_comm_col_size( templ, &nprows );

  mycol_first_rank = me-me%nprows;

  displs[0] = 0;
  for ( i=0; i<nprows; i++ ) {
    PLA_Temp_vector_distant_length(templ,
                              mycol_first_rank+i,
                              length_mv,
                              align_row_mv,
                              &counts[i]);
    counts[i] *= width_mv;
    if ( i>0 ) displs[i] = displs[i-1] + counts[i-1];
  }

  return PLA_SUCCESS;
}

