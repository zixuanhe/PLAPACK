/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/***************************************************************************/

int PLA_Temp_compute_owners_and_sizes_enter(
       int dir,                       /* Direction in which to compute */
       PLA_Template templ,            /* Template for distribution */   
       int global_align,              /* Alignment to template wrt direction */
       int global_length,             /* Length in direction */ 
       int *local_length,             /* Length assigned to this node */
       int *size_first,               /* Size of first block */
       int *owner_first,              /* Owner of first block */
       int *size_last,                /* Size of last block */
       int *owner_last )              /* Owner of last block */

/*----------------------------------------------------------------------------

Purpose : Compute local information for object

IN     dir               Direction in which to compute
IN     templ             Template for distribution
IN     global_align      Alignment to template wrt direction
IN     global_length     Length in direction
OUT    local_length      Length assigned to this node
OUT    size_first        Size of first block
OUT    owner_first       Owner of first block
OUT    size_last         Size of last block
OUT    owner_last        Owner of last block

----------------------------------------------------------------------------*/
{
  return PLA_SUCCESS;
}


int PLA_Temp_compute_owners_and_sizes_exit(
       int dir,                       /* Direction in which to compute */
       PLA_Template templ,            /* Template for distribution */   
       int global_align,              /* Alignment to template wrt direction */
       int global_length,             /* Length in direction */ 
       int *local_length,             /* Length assigned to this node */
       int *size_first,               /* Size of first block */
       int *owner_first,              /* Owner of first block */
       int *size_last,                /* Size of last block */
       int *owner_last )              /* Owner of last block */

/*----------------------------------------------------------------------------

Purpose : Compute local information for object

IN     dir               Direction in which to compute
IN     templ             Template for distribution
IN     global_align      Alignment to template wrt direction
IN     global_length     Length in direction
OUT    local_length      Length assigned to this node
OUT    size_first        Size of first block
OUT    owner_first       Owner of first block
OUT    size_last         Size of last block
OUT    owner_last        Owner of last block

----------------------------------------------------------------------------*/
{
  return PLA_SUCCESS;
}


