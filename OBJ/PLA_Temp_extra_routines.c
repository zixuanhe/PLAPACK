/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/***************************************************************************/

int PLA_Temp_compute_owners_and_sizes(
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
  int
    value = PLA_SUCCESS,
    me, rel_me,
    nprocs, nprows, 
    nb,
    nb_last,
    nb_first, 
    num_blocks_before,
    num_blocks_left,
    cur_owner_first;

  if ( value == PLA_SUCCESS ){
    if ( PLA_COMM_ALL == dir ){
      PLA_Temp_comm_all_rank( templ, &me );
      PLA_Temp_comm_all_size( templ, &nprocs );
      PLA_Temp_nb           ( templ, &nb );
    }
    else if ( PLA_COMM_COL == dir ){
      PLA_Temp_comm_col_rank( templ, &me );
      PLA_Temp_comm_col_size( templ, &nprocs );
      PLA_Temp_nb           ( templ, &nb );
    }
    else if ( PLA_COMM_ROW == dir ){
      PLA_Temp_comm_row_rank( templ, &me );
      PLA_Temp_comm_row_size( templ, &nprocs );
      PLA_Temp_nb           ( templ, &nb );

      PLA_Temp_comm_col_size( templ, &nprows );
    }
    else
      PLA_Abort( "Illegal direction", __LINE__, __FILE__ );

/*  printf("length %d,\n", global_length );
  printf("align  %d,\n", global_align ); */

    *size_first = min( nb - ( global_align%nb ), global_length );
    *size_last = (global_align+global_length) % nb;
    *size_last = min ( ( *size_last == 0 ? nb: *size_last ), global_length );

    if ( PLA_COMM_ROW == dir ) nb *= nprows;

    nb_first = min( nb - ( global_align%nb ), global_length );

    num_blocks_before = global_align / nb;
    *owner_first = num_blocks_before % nprocs;

    global_length -= nb_first;
    *local_length = ( me == *owner_first ? nb_first : 0 );
      
    if ( global_length == 0 ){
      *owner_last = *owner_first;
/*      *size_last  = *size_first; */
    }
    else {
      global_align += nb_first;
      
      cur_owner_first = ( *owner_first+1 )%nprocs;
      rel_me = ( me + nprocs  - cur_owner_first ) % nprocs;

      nb_last = global_length%nb;
      if ( nb_last == 0 ) nb_last = min( global_length, nb );
      global_length -= nb_last;
      
      num_blocks_left = global_length / nb;
      *local_length += ( num_blocks_left / nprocs ) * nb;
      num_blocks_left = num_blocks_left % nprocs;
      
      *owner_last = ( cur_owner_first + num_blocks_left ) % nprocs;
      if ( me == *owner_last ) *local_length += nb_last;
      if ( rel_me < num_blocks_left ) *local_length += nb;
    }
  }
/*  printf("local size %d\n", *local_length );
  printf("size  first %d\n", *size_first );
  printf("owner first %d\n", *owner_first );
  printf("size last %d\n", *size_last );
  printf("owner last %d\n", *owner_last ); */


  return value;
}

/***************************************************************************/

int PLA_Temp_compute_local_length(
       int dir,                       /* Direction in which to compute */
       PLA_Template templ,            /* Template for distribution */   
       int global_align,              /* Alignment to template wrt direction */
       int global_length,             /* Length in direction */ 
       int *local_length )            /* Length assigned to this node */

/*----------------------------------------------------------------------------

Purpose : Compute local length for object

IN     dir               Direction in which to compute
IN     templ             Template for distribution
IN     global_align      Alignment to template wrt direction
IN     global_length     Length in direction
OUT    local_length      Length assigned to this node

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    me, rel_me,
    nprocs, nprows, 
    nb,
    nb_last,
    nb_first, 
    num_blocks_before,
    num_blocks_left,
    cur_owner_first,
    owner_first, owner_last,
    size_first, size_last;

  if ( value == PLA_SUCCESS ){
    if ( PLA_COMM_ALL == dir ){
      PLA_Temp_comm_all_rank( templ, &me );
      PLA_Temp_comm_all_size( templ, &nprocs );
      PLA_Temp_nb           ( templ, &nb );
    }
    else if ( PLA_COMM_COL == dir ){
      PLA_Temp_comm_col_rank( templ, &me );
      PLA_Temp_comm_col_size( templ, &nprocs );
      PLA_Temp_nb           ( templ, &nb );
    }
    else if ( PLA_COMM_ROW == dir ){
      PLA_Temp_comm_row_rank( templ, &me );
      PLA_Temp_comm_row_size( templ, &nprocs );
      PLA_Temp_nb           ( templ, &nb );

      PLA_Temp_comm_col_size( templ, &nprows );
    }
    else
      PLA_Abort( "Illegal direction", __LINE__, __FILE__ );

/*  printf("length %d,\n", global_length );
  printf("align  %d,\n", global_align ); */

    size_first = min( nb - ( global_align%nb ), global_length );
    size_last = (global_align+global_length) % nb;
    size_last = min ( ( size_last == 0 ? nb: size_last ), global_length );

    if ( PLA_COMM_ROW == dir ) nb *= nprows;

    nb_first = min( nb - ( global_align%nb ), global_length );

    num_blocks_before = global_align / nb;
    owner_first = num_blocks_before % nprocs;

    global_length -= nb_first;
    *local_length = ( me == owner_first ? nb_first : 0 );
      
    if ( global_length == 0 ){
      owner_last = owner_first;
    }
    else {
      global_align += nb_first;
      
      cur_owner_first = ( owner_first+1 )%nprocs;
      rel_me = ( me + nprocs  - cur_owner_first ) % nprocs;

      nb_last = global_length%nb;
      if ( nb_last == 0 ) nb_last = min( global_length, nb );
      global_length -= nb_last;
      
      num_blocks_left = global_length / nb;
      *local_length += ( num_blocks_left / nprocs ) * nb;
      num_blocks_left = num_blocks_left % nprocs;
      
      owner_last = ( cur_owner_first + num_blocks_left ) % nprocs;
      if ( me == owner_last ) *local_length += nb_last;
      if ( rel_me < num_blocks_left ) *local_length += nb;
    }
  }

  return value;
}

