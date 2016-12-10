/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Temp_vector_distant_length( PLA_Template templ,
			      int proc,
			      int global_length, 
			      int global_align,
			      int *local_m_proc )

{
  int 
    nb_distr,
    nprocs;

  PLA_Temp_nb( templ, &nb_distr );
  PLA_Temp_comm_all_size( templ, &nprocs );
  
  *local_m_proc = PLA_Local_part( proc, nb_distr, nprocs, 
				   global_length+global_align ) -
                  PLA_Local_part( proc, nb_distr, nprocs, 
				   global_align );

  return PLA_SUCCESS;
} 


int PLA_Local_part( int proc, int nb_distr, int nprocs, int length ) 
{
  int
    value,
    nblocks;

  nblocks = length/nb_distr;

  value   = nblocks/nprocs * nb_distr;
  nblocks = nblocks%nprocs;    /* number of blocks left */

  if ( proc < nblocks )  value += nb_distr;
  if ( proc == nblocks ) value += length%nb_distr;

  return value;
}
