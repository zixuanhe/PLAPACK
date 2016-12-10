/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

extern double sqrt();

double pla2_zabs( PLA_DOUBLE_COMPLEX x )
{
  return( dabs(x.real) + dabs(x.imaginary) );
}

float pla2_cabs ( PLA_COMPLEX x)
{
  return (float)  ( dabs(x.real)+ dabs(x.imaginary) );
}

/****************************************************************************
int PLA_Iamax( PLA_Obj x, PLA_Obj k, PLA_Obj xmax)

  Determines the maximum value (xmax) and the global offset from the first
  element (k) of the distributed object x.  The process is :

  1) get the local result via a call to PLA_Local_iamax().
  2) get the maximum value and rank of the owner thereof via MPI_Allreduce().
  3) the owner processor broadcasts its local index
  4) all compute the global index via a call to
     PLA_Index_local_to_global()
*****************************************************************************/
int PLA_Iamax( PLA_Obj x, PLA_Obj k, PLA_Obj xmax)
{
  int
    *k_buf = NULL;

  MPI_Comm
    comm_all;

  MPI_Datatype
    datatype;

  PLA_Template
    template = NULL;

  int 
    global_alignment,
    local_index,
    myrank,
    nb,
    np,
    zero_or_one,
    length_x, width_x,
    index, index_temp,
    i, i_one = 1,
    value = PLA_SUCCESS;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Iamax_enter( x, k, xmax );

  PLA_Obj_datatype( x, &datatype);
  PLA_Obj_template( x, &template);
  PLA_Obj_global_align (x, &global_alignment);
  PLA_Temp_nb( template, &nb);
  PLA_Temp_zero_or_one ( template, &zero_or_one); 

  PLA_Temp_comm_all_info(template, &comm_all, &myrank, &np);

  if ( datatype == MPI_DOUBLE ){
    struct d_pair{
      double value;
      int    index;
    } my_pair, *all_pairs;
    double *buf_x, *xmax_buf, xmax_temp;

    PLA_Obj_local_buffer( x, (void **) &buf_x );
    PLA_Obj_local_length( x, &length_x );
    PLA_Obj_local_width( x,  &width_x );

    if ( length_x * width_x > 0 ){
      index = PLA_idamax( &length_x, buf_x, &i_one ) - 1;
      my_pair.value = *( buf_x+index );
      my_pair.index = 
	PLA_Index_local_to_global (global_alignment,
				    index,
				    myrank,
				    np,
				    nb,
				    zero_or_one);
    }
    else{
      my_pair.value = 0.0;
      my_pair.index = 0;
    }

    all_pairs = ( struct d_pair * ) PLA_malloc( (size_t) np * sizeof( struct d_pair ) );

    MPI_Allgather( &my_pair, sizeof( struct d_pair ), MPI_CHAR, 
		   all_pairs, sizeof( struct d_pair), MPI_CHAR, comm_all );

    xmax_temp = all_pairs[0].value;
    index_temp = all_pairs[0].index;
    for ( i=1; i<np; i++ ){
      if ( dabs( xmax_temp ) < dabs( all_pairs[i].value ) || 
	   ( dabs( xmax_temp ) == dabs( all_pairs[i].value ) &&
	     index_temp > all_pairs[i].index ) ){
	xmax_temp = all_pairs[i].value;
	index_temp = all_pairs[i].index;
      }
    }

    PLA_Obj_local_buffer( xmax, (void **) &xmax_buf);
    PLA_Obj_local_buffer( k, (void **) &k_buf);

    *xmax_buf = xmax_temp;
    *k_buf = index_temp;

    PLA_free( all_pairs );
  }
  else if ( datatype == MPI_FLOAT ){
    struct d_pair{
      float value;
      int    index;
    } my_pair, *all_pairs;
    float *buf_x, *xmax_buf, xmax_temp;

    PLA_Obj_local_buffer( x, (void **) &buf_x );
    PLA_Obj_local_length( x, &length_x );
    PLA_Obj_local_width( x,  &width_x );

    if ( length_x * width_x > 0 ){
      index = PLA_isamax( &length_x, buf_x, &i_one ) - 1;
      my_pair.value = *( buf_x+index );
      my_pair.index = 
	PLA_Index_local_to_global (global_alignment,
				    index,
				    myrank,
				    np,
				    nb,
				    zero_or_one);
    }
    else{
      my_pair.value = 0.0;
      my_pair.index = 0;
    }

    all_pairs = ( struct d_pair * ) PLA_malloc( (size_t) np * sizeof( struct d_pair ) );

    MPI_Allgather( &my_pair, sizeof( struct d_pair ), MPI_CHAR, 
		   all_pairs, sizeof( struct d_pair), MPI_CHAR, comm_all );

    xmax_temp = all_pairs[0].value;
    index_temp = all_pairs[0].index;
    for ( i=1; i<np; i++ ){
      if ( dabs( xmax_temp ) < dabs( all_pairs[i].value ) || 
	   ( dabs( xmax_temp ) == dabs( all_pairs[i].value ) &&
	     index_temp > all_pairs[i].index ) ){
	xmax_temp = all_pairs[i].value;
	index_temp = all_pairs[i].index;
      }
    }

    PLA_Obj_local_buffer( xmax, (void **) &xmax_buf);
    PLA_Obj_local_buffer( k, (void **) &k_buf);

    *xmax_buf = xmax_temp;
    *k_buf = index_temp;

    PLA_free( all_pairs );
  }
  else if ( datatype == MPI_DOUBLE_COMPLEX ){
    struct d_pair{
      PLA_DOUBLE_COMPLEX value;
      int    index;
    } my_pair, *all_pairs;
    double xmax_norm, xmax_norm_temp;
    PLA_DOUBLE_COMPLEX *buf_x, *xmax_buf, xmax_temp;

    PLA_Obj_local_buffer( x, (void **) &buf_x );
    PLA_Obj_local_length( x, &length_x );
    PLA_Obj_local_width( x,  &width_x );

    if ( length_x * width_x > 0 ){
      index = PLA_izamax( &length_x, buf_x, &i_one ) - 1;
      my_pair.value = *( buf_x+index);
      my_pair.index = 
	PLA_Index_local_to_global (global_alignment,
				    index,
				    myrank,
				    np,
				    nb,
				    zero_or_one);
    }
    else{
      my_pair.value.real = my_pair.value.imaginary = 0.0;
      my_pair.index = 0;
    }

    all_pairs = ( struct d_pair * ) PLA_malloc( (size_t) np * sizeof( struct d_pair ) );

    MPI_Allgather( &my_pair, sizeof( struct d_pair ), MPI_CHAR, 
		   all_pairs, sizeof( struct d_pair), MPI_CHAR, comm_all );

    xmax_temp   = all_pairs[0].value;
    xmax_norm_temp   = pla2_zabs( xmax_temp );
    index_temp = all_pairs[0].index;
    for ( i=1; i<np; i++ ){
      xmax_norm = pla2_zabs( all_pairs[i].value );
      if ( ( xmax_norm_temp < xmax_norm  ) || 
	   ( xmax_norm_temp == xmax_norm &&
	     index_temp > all_pairs[i].index ) ){
	xmax_temp = all_pairs[i].value;
	index_temp = all_pairs[i].index;
	xmax_norm_temp = xmax_norm;
      }
    }

    PLA_Obj_local_buffer( xmax, (void **) &xmax_buf);
    PLA_Obj_local_buffer( k, (void **) &k_buf);

    *xmax_buf = xmax_temp;
    *k_buf = index_temp;

    PLA_free( all_pairs );
  }
  else if ( datatype == MPI_COMPLEX ){
    struct d_pair{
      PLA_COMPLEX value;
      int    index;
    } my_pair, *all_pairs;
    PLA_COMPLEX *buf_x, *xmax_buf, xmax_temp;

    PLA_Obj_local_buffer( x, (void **) &buf_x );
    PLA_Obj_local_length( x, &length_x );
    PLA_Obj_local_width( x,  &width_x );

    if ( length_x * width_x > 0 ){
      index = PLA_icamax( &length_x, buf_x, &i_one ) - 1;
      my_pair.value = *( buf_x+index);
      my_pair.index = 
	PLA_Index_local_to_global (global_alignment,
				    index,
				    myrank,
				    np,
				    nb,
				    zero_or_one);
    }
    else{
      my_pair.value.real = my_pair.value.imaginary = 0.0;
      my_pair.index = 0;
    }

    all_pairs = ( struct d_pair * ) PLA_malloc( (size_t) np * sizeof( struct d_pair ) );

    MPI_Allgather( &my_pair, sizeof( struct d_pair ), MPI_CHAR, 
		   all_pairs, sizeof( struct d_pair), MPI_CHAR, comm_all );

    xmax_temp   = all_pairs[0].value;
    index_temp = all_pairs[0].index;
    for ( i=1; i<np; i++ ){
      if ( pla2_cabs( xmax_temp ) < pla2_cabs( all_pairs[i].value ) || 
	   ( pla2_cabs( xmax_temp ) == pla2_cabs( all_pairs[i].value ) &&
	     index_temp > all_pairs[i].index ) ){
	xmax_temp = all_pairs[i].value;
	index_temp = all_pairs[i].index;
      }
    }

    PLA_Obj_local_buffer( xmax, (void **) &xmax_buf);
    PLA_Obj_local_buffer( k, (void **) &k_buf);

    *xmax_buf = xmax_temp;
    *k_buf = index_temp;

    PLA_free( all_pairs );
  }
  else
    PLA_Abort( "datatype not yet implemented", __LINE__, __FILE__ );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Iamax_exit( x, k, xmax );

  return value;
}


int PLA_Index_local_to_global(int global_alignment,
			      int local_index, 
			      int processor_rank,
			      int number_of_processors,
			      int nb,
			      int zero_or_one)
/*  Returns the global alignment index of the (vector) element
    located at local index "local_index" on processor "processor_rank"
    out of a total of "number_of_processors".  Also needs the
    distribution blocksize "nb" and the global alignment index of
    the object under consideration, and the template's zero_or_one value
    (which specifies whether indices are counted from zero or one.)      */
{
  int
    local_block_index = 0,
    local_block_number =0,
    local_block_size = 0,
    global_block_number = 0,
    owner_first_block = 0,
    first_block_complement_size = 0,
    first_block_corrector = 0,
    global_offset_from_first_element = 0,
    global_offset = 0;

  /* transform the local index and global_alignment so that it is as if 
     zero_or_one equaled zero */
  local_index = local_index - zero_or_one;
  global_alignment = global_alignment - zero_or_one ;

  /* who owns the first block */
  owner_first_block = (global_alignment / nb) % number_of_processors;

  /* Set up the correction for a first partial block.
     First compute the size of the "space left" in the very
     first block due to alignment.  Then the corrector is
     this size (but only for the owner of the first block --
     everyone else's blocks start on nb-boundaries) */
  first_block_complement_size = global_alignment % nb; 
  
  first_block_corrector = (processor_rank == owner_first_block ? 
			   first_block_complement_size : 0);

  /* which of the local blocks (including partial blocks) is this index in? */
  local_block_number = ( local_index + first_block_corrector) / nb;

  /* size of the local block (note correction for indexing from zero )*/
  /* Note that this is the size of the local block that owns the current 
     index, up to and including that index */
  local_block_size = 1 + ((local_index + first_block_corrector) % nb);
  
  /* how does that local block number translate to global block number ? */
  global_block_number = ((local_block_number + (processor_rank<owner_first_block?1:0)) * 
			 number_of_processors ) +
                        (processor_rank - owner_first_block) ;

  /* how far in global index space is this element from the first element ?  

     Rationale : count all full blocks (global_block_num -1)*nb,
     add in the first partial block (nb - first_block_comlement_size),
     then finally add in the size of the (partial) block that holds
     the index in question (local_block_size).  Also note the (-1) 
     to correct for indexing from zero */
  global_offset_from_first_element = global_block_number*nb + local_block_size - 1 
                                       - first_block_complement_size;

  /* transform back so that zero_or_one is taken into account */
  return global_offset_from_first_element + zero_or_one ;
}  














