/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

void pla_array_copy_new( MPI_Datatype, int, int, char *, int, char *, int );

/*----------------------------------------------------------------------*/

int PLA_Pack_from_mv_to_pmv_onto_row( PLA_Obj from, char *sendbuffer )

/************************************************************************
 
  Purpose:  Pack the contents of multivector "from" before calling 
  MPI_Gatherv (in case where target is projected multivector) or
  MPI_Allgatherv (in case where target is duplicated proj. multivector)
  within columns.  
 
   from       -- Object from which to pack.  This object references
                 a multivector

   sendbuffer -- Buffer into which to pack.

   Note: This implementation assumes that the distribution block size 
         is uniform.

   Assumptions: from and to are properly aligned.

*************************************************************************/

{
  int ldim_buf, dummy;

  PLA_Obj_global_width( from, &ldim_buf );

  PLA_Obj_get_local_contents( from, PLA_TRANSPOSE, &dummy, &dummy,
 			      sendbuffer, ldim_buf, 1 );
}

/*----------------------------------------------------------------------*/

int PLA_Unpack_from_mv_to_pmv_onto_row( 
      PLA_Obj from,                         /* original multivector   */
      char *recvbuffer,                     /* recvbuffer[]           */
      int *recvcounts,                      /* recvcounts[]           */
      int *recvdispls,                      /* recvdispls[]
                                               us used by MPI_gatherv */
      PLA_Obj to )                          /* target projected mvec  */

/************************************************************************
 
  Purpose:  Unpack the contents of multivector "from" after calling 
  MPI_Gatherv (in case where target is projected multivector) or
  MPI_Allgatherv (in case where target is duplicated proj. multivector)
  within columns.  
 
   from       -- Object from which to pack.  This object references
                 a multivector

   recvbuffer -- Buffer from which to unpack

   revcounts  -- Sizes of the buffers from the individual nodes

   recvdispls -- Displacements of the buffers from the indiv. nodes

   to         -- object into which to unpack

   Note: This implementation assumes that the distribution block size 
         is uniform.

   Assumptions: from and to are properly aligned.

   Note: data was transposed before calling the gather.

*************************************************************************/

{
  int 
    i,
    nprows, npcols,
    mycol,
    owner, owner_row,
    size_top,
    local_width_to, 
    global_length,  
    local_ldim_to,
    nb_distr,
    typesize;
  char 
    **bufp,
    *buffer_to;     
  PLA_Obj 
    to_cur = NULL, 
    from_cur = NULL;
  PLA_Template 
    templ = NULL;
  MPI_Datatype 
    datatype;

  /* get mesh information */
  PLA_Obj_template( from, &templ );
  PLA_Temp_comm_col_size( templ, &nprows );
  PLA_Temp_comm_row_size( templ, &npcols );
  PLA_Temp_comm_row_rank( templ, &mycol );

  /* get distribution block size */
  PLA_Temp_nb( templ, &nb_distr );
  
  /* length and width of local proj. mvec */
  PLA_Obj_local_width( to, &local_width_to );
  PLA_Obj_global_length( to, &global_length );

  /* datatype of objects */
  PLA_Obj_datatype( from, &datatype );
  MPI_Type_size( datatype, &typesize );

  /* address and local leading dimension of where to put data */
  PLA_Obj_local_buffer( to, (void **) &buffer_to );
  PLA_Obj_local_ldim( to, &local_ldim_to );

  PLA_Obj_view_all( from, &from_cur );

  /* compute pointers to where the buffers from different nodes
     start */
  bufp = (char **) PLA_calloc( nprows, sizeof( char * ) );
  for ( i = 0; i<nprows; i++ ){
    bufp[ i ] = (char *) recvbuffer + recvdispls[ i ] * typesize;
  }

  PLA_Obj_view_all( to, &to_cur );
  PLA_Obj_view_all( from, &from_cur );

  /* Update to_cur and from_cur until end of objects or
     this node owns top */
  while( TRUE ){
    PLA_Obj_split_size( from_cur, PLA_SIDE_TOP, &size_top, &owner );
    if ( 0 == size_top ) break;
    if ( mycol == owner/nprows ) break;
    PLA_Obj_vert_split_2( to_cur, size_top, PLA_DUMMY, &to_cur );
    PLA_Obj_horz_split_2( from_cur, size_top, PLA_DUMMY,
			                      &from_cur );
  }

  if ( size_top != 0 ) {
    /* Copy top to right place.  After this, everything
       is aligned on a block boundary */
    PLA_Obj_local_width( to_cur, &local_width_to );
    owner_row = owner%nprows;
    PLA_Obj_local_buffer( to_cur, (void **)  &buffer_to );
    pla_array_copy_new( datatype, global_length, size_top, 
		   bufp[owner_row], global_length,
		   buffer_to, local_ldim_to );
    bufp[ owner_row ] += size_top * typesize * global_length;
/*    buffer_to += size_top * typesize * global_length; */
    buffer_to += size_top * typesize * local_ldim_to;

    /* interleave rest of blocks */
    local_width_to -= size_top;
    for ( ;local_width_to>0; local_width_to-=nb_distr ){
      owner_row = (owner_row+1)%nprows;
      size_top = ( local_width_to < nb_distr ? local_width_to : nb_distr );
      pla_array_copy_new( datatype, global_length, size_top,
		      bufp[ owner_row ], global_length,
	              buffer_to, local_ldim_to );
      bufp[ owner_row ] += size_top * typesize * global_length;
      buffer_to += size_top * typesize * local_ldim_to;
    }
  }

  PLA_free( bufp );

  PLA_Obj_free( &to_cur );
  PLA_Obj_free( &from_cur );

  return ( PLA_SUCCESS );
}  

/*----------------------------------------------------------------------*/

int PLA_Unpack_from_mv_to_pmv_onto_col( 
      PLA_Obj from, 
      char *recvbuffer, 
      int *recvcounts, 
      int *recvdispls, 
      PLA_Obj to )

/************************************************************************
 
  Purpose:  Unpack the contents of multivector "from" after calling 
  MPI_Gatherv (in case where target is projected multivector) or
  MPI_Allgatherv (in case where target is duplicated proj. multivector)
  within rows.
 
   from       -- Object from which to pack.  This object references
                 a multivector

   recvbuffer -- Buffer from which to unpack

   revcounts  -- Sizes of the buffers from the individual nodes

   recvdispls -- Displacements of the buffers from the indiv. nodes

   to         -- object into which to unpack

   Note: This implementation assumes that the distribution block size 
         is uniform.

   Assumptions: from and to are properly aligned.

   Note: data was transposed before calling the gather.

*************************************************************************/

{
  int 
    i,
    nprows, npcols,
    myrow,
    owner, owner_row, owner_col,
    size_top,
    local_length_to, global_width,
    local_ldim_to,
    nb_distr,
    typesize,
    *ldims;
  char
    **bufp,
    *buffer_to;     
  PLA_Obj 
    to_cur = NULL, 
    from_cur = NULL;
  PLA_Template 
    templ = NULL;
  MPI_Datatype 
    datatype;

  /* get mesh information */
  PLA_Obj_template( from, &templ );
  PLA_Temp_comm_col_size( templ, &nprows );
  PLA_Temp_comm_row_size( templ, &npcols );
  PLA_Temp_comm_col_rank( templ, &myrow );

  /* get distribution block size */
  PLA_Temp_nb( templ, &nb_distr );
  
  /* length and width of local proj. mvec */
  PLA_Obj_local_length( to, &local_length_to );
  PLA_Obj_global_width( to, &global_width );

  /* datatype of objects */
  PLA_Obj_datatype( from, &datatype );
  MPI_Type_size( datatype, &typesize );

  /* address and local leading dimension of where to put data */
  PLA_Obj_local_buffer( to, (void **) &buffer_to );
  PLA_Obj_local_ldim( to, &local_ldim_to );

  /* compute leading dimensions of the recv buffers */
  ldims = (int *) PLA_calloc( npcols, sizeof( int ) );
  for (i=0; i<npcols; i++) ldims[i] = recvcounts[i]/global_width;

  PLA_Obj_view_all( from, &from_cur );

  /* compute pointers to where the buffers from different nodes
     start */
  bufp = (char **) PLA_calloc( npcols, sizeof( char * ) );
  for ( i = 0; i<npcols; i++ ){
    bufp[ i ] = (char *) recvbuffer + recvdispls[ i ] * typesize;
  }

  PLA_Obj_view_all( to, &to_cur );
  PLA_Obj_view_all( from, &from_cur );

  /* Update to_cur and from_cur until end of objects or
     this node owns top */
  while( TRUE ){
    PLA_Obj_split_size( from_cur, PLA_SIDE_TOP, &size_top, &owner );
    if ( 0 == size_top ) break;
    if ( myrow == ( owner_row = owner%nprows ) ) break;
    PLA_Obj_horz_split_2( to_cur, size_top, PLA_DUMMY,
			                   &to_cur );
    PLA_Obj_horz_split_2( from_cur, size_top, PLA_DUMMY,
			                      &from_cur );
  }

  if ( size_top != 0 ) {
    /* Copy top to right place.  After this, everything
       is aligned on a block boundary */
    PLA_Obj_local_length( to_cur, &local_length_to );
    owner_col = owner/nprows;
    PLA_Obj_local_buffer( to_cur, (void **)  &buffer_to );
    pla_array_copy_new( datatype, size_top, global_width, 
		   bufp[owner_col], ldims[ owner_col ],
		   buffer_to, local_ldim_to );
    bufp[ owner_col ] += size_top * typesize;
    buffer_to += size_top * typesize;

    /* interleave rest of blocks */
    local_length_to -= size_top;
    for ( ;local_length_to>0; local_length_to-=nb_distr ){
      owner_col = (owner_col+1)%npcols;
      size_top = ( local_length_to < nb_distr ? local_length_to : nb_distr );
      pla_array_copy_new( datatype, size_top, global_width, 
		      bufp[ owner_col ], ldims[ owner_col ],
	              buffer_to, local_ldim_to );
      bufp[ owner_col ] += size_top * typesize;
      buffer_to += size_top * typesize;
    }
  }

  PLA_free( bufp );
  PLA_free( ldims );

  PLA_Obj_free( &to_cur );
  PLA_Obj_free( &from_cur );

  return ( PLA_SUCCESS );
}  
                  
/*----------------------------------------------------------------------*/

int PLA_Pack_from_pmv_onto_col_to_mv( 
      PLA_Obj from,              
      char *sendbuffer,          
      int *sendcounts,            
      int *senddispls, 
      PLA_Obj to )

/************************************************************************

   Purpose: pack data for scatter.
   Assumptions: from and to are properly aligned.

   from       -- Object from which to pack.  This object references
                 a projected multivector or matrix that exists entirely 
                 within one column of nodes, and is only called by nodes
                 that own parts of these columns.

   sendbuffer -- Buffer into which to pack.

   sendcounts -- Array of lengths to be passed to MPI_Scatterv.

   senddispls -- Array of displacements to be passed to MPI_Scatterv.

   to         -- Target object.

   Note: This implementation assumes that the distribution block size 
         is uniform.

************************************************************************/

{
  int 
    i,              niters,
    nprows,         npcols,
    myrow,
    owner,          owner_col,
    size_top,
    *ldims,
    local_length_from,
    global_width,
    local_ldim_from,
    nb_distr,
    typesize;
  char 
    **bufp,
    *buffer_from;     
  PLA_Obj 
    to_cur = NULL, 
    from_cur = NULL;
  PLA_Template 
    templ = NULL;
  MPI_Datatype 
    datatype;

  /* Extract relevant information from the template */
  PLA_Obj_template( from, &templ );
  PLA_Temp_comm_col_rank( templ, &myrow );
  PLA_Temp_comm_col_size( templ, &nprows );
  PLA_Temp_comm_row_size( templ, &npcols );
  PLA_Temp_nb( templ, &nb_distr );
  
  /* Extract relevant information from the source object */
  PLA_Obj_datatype    ( from, &datatype );
  PLA_Obj_global_width( from, &global_width );
  PLA_Obj_local_ldim  ( from, &local_ldim_from );

  /* Determine size of the data items */
  MPI_Type_size( datatype, &typesize );

  /* Create an array in which to store the leading dimensions
     of each of the subarrays of the send buffer. */

  ldims = (int *) PLA_malloc( (size_t) npcols * sizeof( int ) );
  for ( i=0; i<npcols; i++ ) ldims[i] = sendcounts[i] / global_width;

  /* Create an array of pointers that will track where to put
     the next contribution to be sent to node ( myrow, i ) */

  bufp = (char **) PLA_malloc( (size_t) npcols * sizeof( char * ) );
  for ( i = 0; i<npcols; i++ )
    bufp[ i ] = (char *) sendbuffer + senddispls[ i ] * typesize;


  /* Determine  to = /   *    \
                     \ to_cur /
     where the first item of to_cur is the first item of "to"
     assigned to this row of nodes. Upon exit of loop, owner
     will equal the node index (relative to all nodes) that owns
     the top of to_cur.  from_cur is taken to view the corresponding
     part of "from" */

  PLA_Obj_view_all( to, &to_cur );
  PLA_Obj_view_all( from, &from_cur );

  while( TRUE ){
    PLA_Obj_split_size( to_cur, PLA_SIDE_TOP, &size_top, &owner );
    if ( 0 == size_top ) break;
    if ( myrow == owner%nprows ) break;
    PLA_Obj_horz_split_2( from_cur, size_top, PLA_DUMMY,
			                      &from_cur );
    PLA_Obj_horz_split_2( to_cur, size_top, PLA_DUMMY,
			                    &to_cur );
  }

  if ( size_top != 0 ) {
    /* The following code is equivalent to
       PLA_Obj_split_size( from_cur, PLA_SIDE_TOP, &size_top, &owner );
       PLA_Obj_split_2( from_cur, size_top,  &from_1,
                                             &from_cur );
       PLA_Obj_split_2( to_cur, size_top,    &to_1,
                                             &to_cur );
       --- pack data in from_1 into appropriate part of send buffer ---

       After this, the buffer_from is updated so that it points to
       a clean subvector of length nb, so that a loop can be used
       to copy the local buffer (pointed to by buffer_from) to the
       subarrays of the sendbuffer. */

    /* local_length_from is the total number of rows to be packed */
    PLA_Obj_local_length( from_cur, &local_length_from);

    /* owner_col equals the column index of the node to which the
       first "size_top" rows must be sent, where size_top was 
       chosen to extend to the next processor boundary */
       
    owner_col = owner/nprows;

    /* Let buffer_from point to the local buffer to be packed */

    PLA_Obj_local_buffer( from_cur, (void **)  &buffer_from );

    /* Pack the first size_top rows to the appropriate subarray of
       the send buffer */

    pla_array_copy_new( datatype, size_top, global_width, 
		   buffer_from, local_ldim_from,
  		   bufp[ owner_col ], ldims[ owner_col ] ); 

    /* Update the pointer into that subarray of the sendbuffer */

    bufp[ owner_col ] += size_top * typesize;

    /* Update the pointer to the buffer from which to copy */

    buffer_from += size_top * typesize;

    /* Update the number of rows yet to be packed */

    local_length_from -= size_top;

    /* Pack the rest of the rows in a loop, taking advantage of the
       fact that we know which node owns the next part of the 
       data to be packed, and that that data starts on a block boundary. */

#ifdef PLA_OMP
    /* niters = ceil( local_length_from / nb_distr ) */
    niters = (local_length_from+nb_distr-1) / nb_distr; 

    /* printf("PLA_Pack_from_pmv_onto_col_to_mv -- new loop with niters = %d\n", niters); */

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
    for (i = 0 ; i < niters; i++ ){
      int my_local_length_from, my_owner_col, my_size_top;
      char *my_buffer_from, *my_bufp;

      my_local_length_from = local_length_from - (i * nb_distr);

      my_owner_col = (owner_col+i+1)%npcols;

      my_size_top = ( my_local_length_from < nb_distr ? my_local_length_from : nb_distr );

      my_buffer_from = buffer_from + (i * nb_distr * typesize);

      my_bufp = bufp[ my_owner_col ] + ((i/npcols) * nb_distr * typesize);

      pla_array_copy_new( datatype, my_size_top, global_width, 
		       my_buffer_from, local_ldim_from,
		       my_bufp, ldims[ my_owner_col ] ); 
    }

    /* printf("PLA_Pack_from_pmv_onto_col_to_mv -- loop complete\n"); */
#else
    for ( ;local_length_from>0; local_length_from-=nb_distr ){
      owner_col = (owner_col+1)%npcols;

      size_top = ( local_length_from < nb_distr ? local_length_from : nb_distr );
      pla_array_copy_new( datatype, size_top, global_width, 
		       buffer_from, local_ldim_from,
		       bufp[ owner_col ], ldims[ owner_col ] ); 

      bufp[ owner_col ] += size_top * typesize;

      buffer_from += size_top * typesize;
    }
#endif
  }

  /* Free the temporary arrays */
  PLA_free( bufp );
  PLA_free( ldims );

  /* Free the temporary objects */
  PLA_Obj_free( &to_cur );
  PLA_Obj_free( &from_cur );

  return PLA_SUCCESS;
}  

                  
/*----------------------------------------------------------------------*/

int PLA_Pack_from_pmv_onto_row_to_mv( 
      PLA_Obj from,                        /* original projected mvec */
      char *sendbuffer,                     /* sendbuffer[]           */
      int *sendcounts,                      /* sendcounts[]           */
      int *senddispls,                      /* senddispls[]
                                               us used by MPI_scatterv */
      PLA_Obj to )                          /* target multivector     */

/************************************************************************
 
  Purpose:  Pack the contents of projected mvector "from" before calling 
  MPI_Scatterv (in case where source is projected multivector) 
  within columns.  
 
   from       -- Object from which to pack.  This object references
                 a projected multivector

   sendbuffer -- Buffer into which to pack

   sendcounts  -- Sizes of the buffers for the individual nodes

   senddispls -- Displacements of the buffers for the indiv. nodes

   to         -- object into which to unpack

   Note: This implementation assumes that the distribution block size 
         is uniform.

   Assumptions: from and to are properly aligned.

   Note: data will be transposed after calling the scatter.

*************************************************************************/

{
  int 
    i, niters,
    nprows, npcols,
    mycol,
    owner, owner_row,
    size_top,
    local_width_from, 
    global_length,  
    local_ldim_from,
    nb_distr,
    typesize;
  char 
    **bufp,
    *buffer_from;     
  PLA_Obj 
    to_cur = NULL, 
    from_cur = NULL;
  PLA_Template 
    templ = NULL;
  MPI_Datatype 
    datatype;

  /* get mesh information */
  PLA_Obj_template( from, &templ );
  PLA_Temp_comm_col_size( templ, &nprows );
  PLA_Temp_comm_row_size( templ, &npcols );
  PLA_Temp_comm_row_rank( templ, &mycol );

  /* get distribution block size */
  PLA_Temp_nb( templ, &nb_distr );
  
  /* length and width of local proj. mvec */
  PLA_Obj_local_width( from, &local_width_from );
  PLA_Obj_global_length( from, &global_length );

  /* datatype of objects */
  PLA_Obj_datatype( from, &datatype );
  MPI_Type_size( datatype, &typesize );

  /* address and local leading dimension of where to find data */
  PLA_Obj_local_buffer( from, (void **) &buffer_from );
  PLA_Obj_local_ldim( from, &local_ldim_from );

  PLA_Obj_view_all( to, &to_cur );

  /* compute pointers to where the buffers from different nodes
     start */
  bufp = (char **) PLA_calloc( nprows, sizeof( char * ) );
  for ( i = 0; i<nprows; i++ ){
    bufp[ i ] = (char *) sendbuffer + senddispls[ i ] * typesize;
  }

  PLA_Obj_view_all( from, &from_cur );
  PLA_Obj_view_all( to, &to_cur );

  /* Update to_cur and from_cur until end of objects or
     this node owns top */
  while( TRUE ){
    PLA_Obj_split_size( to_cur, PLA_SIDE_TOP, &size_top, &owner );
    if ( 0 == size_top ) break;
    if ( mycol == owner/nprows ) break;
    PLA_Obj_vert_split_2( from_cur, size_top, PLA_DUMMY, &from_cur );
    PLA_Obj_horz_split_2( to_cur,   size_top, PLA_DUMMY,
			                      &to_cur );
  }

  if ( size_top != 0 ) {
    /* Copy top to right place.  After this, everything
       is aligned on a block boundary */
    PLA_Obj_local_width( from_cur, &local_width_from );
    owner_row = owner%nprows;
    PLA_Obj_local_buffer( from_cur, (void **)  &buffer_from );
    pla_array_copy_new( datatype, global_length, size_top, 
		   buffer_from, local_ldim_from,
		   bufp[owner_row], global_length );
    bufp[ owner_row ] += size_top * typesize * global_length;
/*    buffer_to += size_top * typesize * global_length; */
    buffer_from += size_top * typesize * local_ldim_from;

    /* interleave rest of blocks */
    local_width_from -= size_top;
#ifdef PLA_OMP
    /* niters = ceil( local_width_from / nb_distr ) */
    niters = (local_width_from+nb_distr-1) / nb_distr; 

    /* printf("PLA_Pack_from_pmv_onto_row_to_mv -- new loop with niters = %d\n", niters); */

#ifdef _OPENMP
#pragma omp parallel for default(shared) schedule(static)
#endif
     for (i = 0 ; i < niters; i++ ){
      int my_local_width_from, my_owner_row, my_size_top;
      char *my_buffer_from, *my_bufp;

      my_local_width_from = local_width_from - (i * nb_distr);

      my_owner_row = (owner_row+i+1)%nprows;

      my_size_top = ( my_local_width_from < nb_distr ? my_local_width_from : nb_distr );

      my_buffer_from = buffer_from + (i * nb_distr * typesize * local_ldim_from);

      my_bufp = bufp[ my_owner_row ] + ((i/nprows) * nb_distr * typesize * global_length);

      pla_array_copy_new( datatype, global_length, my_size_top, 
		       my_buffer_from, local_ldim_from,
		       my_bufp , global_length ); 
    }

    /* printf("PLA_Pack_from_pmv_onto_row_to_mv -- loop complete\n"); */
#else
    for ( ;local_width_from>0; local_width_from-=nb_distr ){
      owner_row = (owner_row+1)%nprows;
      size_top = ( local_width_from < nb_distr ? local_width_from : nb_distr );
      pla_array_copy_new( datatype, global_length, size_top,
	              buffer_from, local_ldim_from,
		      bufp[ owner_row ], global_length );
      bufp[ owner_row ] += size_top * typesize * global_length;
      buffer_from += size_top * typesize * local_ldim_from;
    }
#endif
  }

  PLA_free( bufp );

  PLA_Obj_free( &to_cur );
  PLA_Obj_free( &from_cur );

  return ( PLA_SUCCESS );
}  
                  

/*----------------------------------------------------------------------*/

int PLA_Unpack_from_pmv_onto_row_to_mv( PLA_Obj to, char *recvbuffer )

/************************************************************************
 
  Purpose:  Unpack the contents of multivector "to" after calling 
  MPI_Scatterv 
  within columns.  
 
   to          -- Object to which to unpack.  This object references
                 a multivector

   recvbuffer -- Buffer from which to unpack.

   Note: This implementation assumes that the distribution block size 
         is uniform.

   Assumptions: from and to are properly aligned.

*************************************************************************/

{
  int ldim_buf, length, width;

  PLA_Obj_global_width( to, &ldim_buf );

  PLA_Obj_local_length( to, &length );
  PLA_Obj_local_width( to, &width );

  PLA_Obj_set_local_contents( PLA_TRANSPOSE, width, length,
 			      recvbuffer, ldim_buf, 1, to );
}
