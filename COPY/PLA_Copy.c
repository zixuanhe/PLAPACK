/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include "PLA_Copy_new.h"

int report_case = FALSE;

/*----------------------------------------------------------------------*/

int PLA_Copy( PLA_Obj Obj_from, PLA_Obj Obj_to )

/************************************************************************
 
  Copy contents of Obj_from to Obj_to

*************************************************************************/

{
  int 
    value,
    objtype_from, proj_onto, 
    owner_col, owner_row,
    ooc_from, ooc_to;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Copy_enter( Obj_from, Obj_to );

  PLA_Obj_ooc( Obj_from, &ooc_from );
  PLA_Obj_ooc( Obj_to, &ooc_to );

  if ( ooc_from || ooc_to ){
    return PLA_OOC_Copy( Obj_from, Obj_to );
  }

  PLA_Obj_objtype( Obj_from, &objtype_from );
  
  switch( objtype_from ){
  case PLA_MSCALAR: 
    PLA_Copy_from_msc( Obj_from, Obj_to );
    break;
/*  case PLA_VECTOR:   */
  case PLA_PMVECTOR: 
    PLA_Obj_project_onto( Obj_from, &proj_onto );
    switch( proj_onto ){
    case PLA_PROJ_ONTO_ROW: 
    case PLA_PROJ_ONTO_COL:
      PLA_Obj_owner_col( Obj_from, &owner_col );
      PLA_Obj_owner_row( Obj_from, &owner_row );
      if ( owner_col == PLA_ALL_COLS ||
		      owner_row == PLA_ALL_ROWS )
	PLA_Copy_from_dpmv( Obj_from, Obj_to );
      else
	PLA_Copy_from_pmv( Obj_from, Obj_to );
      break;
    }
    break;
  case PLA_MVECTOR: 
    PLA_Copy_from_mv( Obj_from, Obj_to );
    break;
  case PLA_MATRIX:  
    PLA_Copy_from_matrix( Obj_from, Obj_to );
    break;
  default:
    PLA_Abort ( "Obj_from in PLA_Copy has unknown objtype", __LINE__, __FILE__ );
    exit( 0 );
  }
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Copy_exit( Obj_from, Obj_to );

  return value;
}

int PLA_Copy_from_matrix( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  int objtype_to, proj_onto, owner_col, owner_row;

  PLA_Obj_objtype( Obj_to, &objtype_to );
  
  switch( objtype_to ){
  case PLA_MSCALAR: 
    PLA_Copy_from_matrix_to_msc( Obj_from, Obj_to );
    break;
/*  case PLA_VECTOR:  */
  case PLA_PMVECTOR: 
    PLA_Obj_project_onto( Obj_to, &proj_onto );
    switch( proj_onto ){
    case PLA_PROJ_ONTO_ROW: 
    case PLA_PROJ_ONTO_COL:
      PLA_Obj_owner_col( Obj_to, &owner_col );
      PLA_Obj_owner_row( Obj_to, &owner_row );
      if ( owner_col == PLA_ALL_COLS ||
	  owner_row == PLA_ALL_ROWS )
	PLA_Copy_from_matrix_to_dpmv( Obj_from, Obj_to );
      else
	PLA_Copy_from_matrix_to_pmv( Obj_from, Obj_to );
      break;
    }
    break;
  case PLA_MVECTOR: 
    PLA_Copy_from_matrix_to_mv( Obj_from, Obj_to );
    break;
  case PLA_MATRIX:  
    PLA_Copy_from_matrix_to_matrix( Obj_from, Obj_to );
    break;
  default:
    printf("Obj_to in PLA_Copy has unknown objtype\n");
    exit( 0 );
  }
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_matrix_to_msc( PLA_Obj Obj_from, PLA_Obj Obj_to )

/************************************************************************
 
  Copy from matrix to multiscalar

*************************************************************************/
{
  int 
    done = FALSE,
    global_length,    global_width,
    size_top,         size_left,
    owner_top,        owner_left;

  PLA_Obj_global_length( Obj_from, &global_length );
  PLA_Obj_split_size( Obj_from, PLA_SIDE_TOP, &size_top, &owner_top );

  PLA_Obj_global_width ( Obj_from, &global_width );
  PLA_Obj_split_size( Obj_from, PLA_SIDE_LEFT, &size_left, &owner_left );
  
  if ( global_length == size_top && global_width == size_left ){
    /* Data exists on one node.  */
    PLA_Obj
      from_cur = NULL;

    PLA_Obj_view_all( Obj_from, &from_cur );
    PLA_Obj_objtype_cast( from_cur, PLA_MSCALAR );
    PLA_Copy_from_msc_to_msc( from_cur, Obj_to );

    PLA_Obj_free( &from_cur );

    done = TRUE;
  }
  else{
    PLA_Obj
      local_contr = NULL;

    PLA_Mscalar_create_conf_to( Obj_from, PLA_ALL_ROWS, PLA_ALL_COLS, 
			        &local_contr );

    PLA_Obj_set_to_zero( local_contr );
    PLA_Copy_local_part_from_matrix_to_msc( Obj_from, local_contr ); 
    PLA_Reduce( local_contr, MPI_SUM, Obj_to );

    PLA_Obj_free( &local_contr );

    done = TRUE;
  }

  if ( !done ){
    if ( report_case ) {
      printf("PLA_Copy_from_matrix_to_msc not yet implemented\n");
      report_case = FALSE;
    }
    return PLA_Copy_old( Obj_from, Obj_to );
  }


  return PLA_SUCCESS;
}


int PLA_Copy_from_matrix_to_matrix( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  PLA_Obj from_cur = NULL,   from_1 = NULL,
          to_cur   = NULL,   to_1   = NULL;
  PLA_Template templ;
  int proj_onto_from,     proj_onto_to,
      align_col_from,     align_col_to,
      align_row_from,     align_row_to,
      size_from,          size_to,         size,
      owner_from,         owner_to,
      local_length,       local_width,
      myrow,              mycol,
      nb_distr,           dummy,
      type_size,          done;
  MPI_Comm comm_row = MPI_COMM_NULL,   comm_col = MPI_COMM_NULL;
  MPI_Datatype datatype;
  MPI_Status status;
  void *buffer_temp;

  PLA_Obj_get_orientation( Obj_from, &proj_onto_from );
  PLA_Obj_get_orientation( Obj_to,   &proj_onto_to );

  done = FALSE;

  PLA_Obj_global_align_col( Obj_from, &align_col_from );
  PLA_Obj_global_align_col( Obj_to, &align_col_to );
    
  PLA_Obj_global_align_row( Obj_from, &align_row_from );
  PLA_Obj_global_align_row( Obj_to, &align_row_to );
  if ( proj_onto_from == proj_onto_to ){
    if ( proj_onto_from == PLA_PROJ_ONTO_ROW ){
      if ( align_col_from == align_col_to ) {
	if ( align_row_from == align_row_to ){ 
	  PLA_Local_copy( Obj_from, Obj_to );
	  done = TRUE;
	}
	else {  /* align_row_from != align_row_to */
	  PLA_Obj_local_width( Obj_from, &local_width );
	  if ( 0 != local_width ) {
	    PLA_Obj_template( Obj_from, &templ );
	    PLA_Temp_nb( templ, &nb_distr );
	    PLA_Temp_comm_col( templ, &comm_col );
	    PLA_Temp_comm_col_rank( templ, &myrow );
	    PLA_Temp_comm_row_rank( templ, &mycol );

	    PLA_Obj_datatype( Obj_from, &datatype );
	    MPI_Type_size( datatype, &type_size );
	    buffer_temp = (void *) PLA_malloc( (size_t) nb_distr * local_width * type_size );
	    PLA_Obj_view_all( Obj_from, &from_cur );
	    PLA_Obj_view_all( Obj_to,   &to_cur );

	    while ( TRUE ){
	      PLA_Obj_split_size( from_cur, PLA_SIDE_TOP, &size_from, &owner_from );
	      PLA_Obj_split_size( to_cur,   PLA_SIDE_TOP, &size_to,   &owner_to);

	      if ( 0 == ( size = min( size_from, size_to ) ) ) break;
	      
	      PLA_Obj_horz_split_2( from_cur, size,   &from_1,
				                      &from_cur );
	      PLA_Obj_horz_split_2( to_cur,   size,   &to_1,
                                                      &to_cur );

	      if ( myrow == owner_from && owner_from == owner_to ){
		PLA_Local_copy( from_1, to_1 );
	      }
	      else{
		if ( myrow == owner_from ) {
		  PLA_Obj_get_local_contents( from_1, PLA_NO_TRANS, 
                                               &dummy, &dummy,
					       buffer_temp, size, 1 );  

		  MPI_Send( BF( buffer_temp ), size * local_width, datatype,
			    owner_to, 0, comm_col );
		}
		if ( myrow == owner_to ) {
		  MPI_Recv( BF( buffer_temp ), size * local_width, datatype,
			    owner_from, MPI_ANY_TAG, comm_col, &status );

		  PLA_Obj_set_local_contents( PLA_NO_TRANS, size, local_width,
					       buffer_temp, size, 1, to_1 ); 
		}
	      }
	    }
	    PLA_free( buffer_temp );
	  }
	  
	  done = TRUE;
	}
      }
    }
    else if ( proj_onto_from == PLA_PROJ_ONTO_COL ){
      if ( align_col_from == align_col_to ) {
	if ( align_row_from == align_row_to ){ 
	  PLA_Local_copy( Obj_from, Obj_to );
	  done = TRUE;
	}
	else {  /* align_row_from != align_row_to */
	}
      }
    }
    else {

      if ( align_row_from == align_row_to ) {

      }
    }
    PLA_Obj_free( &from_cur );
    PLA_Obj_free( &from_1 );
    PLA_Obj_free( &to_cur );
    PLA_Obj_free( &to_1 );

  }
  else {
    int 
      size, proj_onto;
    PLA_Obj
      temp = NULL;

    PLA_Obj_project_onto( Obj_from, &proj_onto );

    if ( proj_onto == PLA_PROJ_ONTO_COL )
      PLA_Obj_global_width( Obj_from, &size );
    else
      PLA_Obj_global_width( Obj_to, &size );

    PLA_Mvector_create_conf_to( Obj_from, size, &temp );

    PLA_Copy_from_matrix_to_mv( Obj_from, temp );
    PLA_Copy_from_mv_to_matrix( temp, Obj_to );

    PLA_Obj_free( &temp );

    done = TRUE; 
  }
      
  if ( !done ){
    if ( report_case ) {
      printf("PLA_Copy_from_matrix_to_matrix not yet implemented\n");
      report_case = FALSE;
    }
    return PLA_Copy_old( Obj_from, Obj_to );
  }

}



/*----------------------------------------------------------------------*/

int PLA_Copy_from_matrix_to_mv( PLA_Obj Obj_from, PLA_Obj Obj_to )

/************************************************************************
 
  Copy from matrix to multivector

*************************************************************************/
{
  int proj_onto, 
      align_row_from, align_col_from,       
      align_row_to, align_col_to,
      me,  myrow, mycol, myrow_first_rank, mycol_first_rank,
      nprows,     npcols,
      size, 
      owner_row, owner_col,
      *ldims,
      *sendcounts,
      *senddispls,
      length_mv,
      local_length_mv,
      i;
  void *sendbuf,
       *recvbuf;
  PLA_Obj   A_cur = NULL,        A_1 = NULL,     A_1_copy = NULL,
            mv_cur = NULL,       mv_1 = NULL,    mv_1_copy = NULL;
  PLA_Template templ = NULL;
  MPI_Comm  comm = MPI_COMM_NULL;
  MPI_Datatype datatype;

  PLA_Obj_project_onto( Obj_from, &proj_onto );
  PLA_Obj_datatype( Obj_from, &datatype );

  switch( proj_onto ){
  case PLA_PROJ_ONTO_COL:
    /* Case: matrix is marked projected onto columns. */
    PLA_Obj_global_align_row( Obj_from, &align_row_from );
    PLA_Obj_global_align_row( Obj_to,   &align_row_to );
    if ( align_row_from == align_row_to ){
      /* Case: matrix and multivector a perfectly aligned. */

      /* Get global and local length of the objects */
      PLA_Obj_global_length( Obj_to, &length_mv );
      PLA_Obj_local_length( Obj_to, &local_length_mv );

      /* Extract me, mycol, nprows, npcols, and row communicator */
      PLA_Obj_template( Obj_from, &templ );
      PLA_Temp_comm_row_rank( templ, &mycol );  /* mycol = col index in mesh */
      PLA_Temp_comm_all_rank( templ, &me );     /* me = index in all nodes   */
      PLA_Temp_comm_row( templ, &comm );        /* row communicator          */
      PLA_Temp_comm_col_size( templ, &nprows ); /* number of rows    in mesh */
      PLA_Temp_comm_row_size( templ, &npcols ); /* number of columns in mesh */

      /* ldims[ i ] will hold the number of rows of the target multivector
         that node ( myrow, i ) holds.                                       */

      ldims      = (int *) PLA_calloc( npcols, sizeof(int) ); 

      /* sendcounts[ i ] will hold the number of bytes to be sent to 
         node ( myrow, i ) ( = ldims[i] * size, where size will be 
         the width of the current panel being redistributed)                 */

      sendcounts = (int *) PLA_calloc( npcols, sizeof(int) );

      /* senddispls[ i ] will hold the displacement where the data to be
         sent to node ( myrow, i ) starts                                    */

      senddispls = (int *) PLA_calloc( npcols, sizeof(int) );

      /* Determine the index (relative to all nodes) of the first
         node in this node's row of nodes                                    */

      myrow_first_rank = me%nprows;

      /* compute ldims[i] = the local length of the target multivector on
         node ( myrow, i )                                                   */

      for ( i=0; i<npcols; i++ )
        PLA_Temp_vector_distant_length(templ,
                                  myrow_first_rank+i*nprows, 
                                  length_mv,
                                  align_row_to, 
                                  &ldims[i]);

      /* Take a views into Obj_from and Obj_to so that we can march
         through them slicing off column panels of Obj_from that all 
         exist within one column of nodes                                    */

      PLA_Obj_view_all( Obj_from, &A_cur );
      PLA_Obj_view_all( Obj_to,   &mv_cur );

      while( TRUE ){
	/* Determine width of column panel that exists all on one column of
           nodes                                                             */

	PLA_Obj_split_size_to_next_proc( A_cur, PLA_SIDE_LEFT, &size, 
					 &owner_col ); 

	if ( 0 == size ) break;

	/* Split off current panel to be redistributed                       */


	PLA_Obj_vert_split_2( A_cur,  size,    &A_1,  &A_cur );
	PLA_Obj_vert_split_2( mv_cur, size,    &mv_1, &mv_cur );

	/* Create a workspace into which to pack the current panel */
        PLA_Matrix_create_conf_to( A_1, &A_1_copy );

	/* Compute the sendcounts and displacements for the scatter 
           operation                                                         */
	senddispls[0] = 0;
	for ( i=0; i<npcols; i++ ) {
	  sendcounts[i] = ldims[i] * size;
	  if ( i>0 ) senddispls[i] = senddispls[i-1] + sendcounts[i-1];
	}

	/* If this node owns part of current panel, pack                     */

	if ( mycol == owner_col ){
	  PLA_Obj_local_buffer( A_1_copy, &sendbuf );

	  PLA_Pack_from_pmv_onto_col_to_mv( 
             A_1, sendbuf, sendcounts, senddispls, mv_1 );
	}

	/* Create a workspace into which to receive this node's 
           part of the multivector */

	PLA_Mvector_create_conf_to( mv_1, size, &mv_1_copy );
	PLA_Obj_local_buffer( mv_1_copy, &recvbuf );

	/* Scatter the data within rows.                                     */

	MPI_Scatterv( BF( sendbuf ), sendcounts, senddispls, datatype, 
                      BF( recvbuf ), local_length_mv * size, datatype, 
                      owner_col, comm );

	/* Copy into the multivector.  (Notice: we cannot receive 
           straight into mv_1 since the local leading dimension
           may not equal the number of rows being received.)                 */
	PLA_Local_copy( mv_1_copy, mv_1 );
      }

      /* Free temporary objects                                              */
      PLA_Obj_free( &A_cur );        
      PLA_Obj_free( &A_1 );
      PLA_Obj_free( &mv_cur );       
      PLA_Obj_free( &mv_1 );
      PLA_Obj_free( &A_1_copy );     
      PLA_Obj_free( &mv_1_copy );

      /* Free temporary arrays                                               */
      PLA_free( ldims );
      PLA_free( sendcounts );
      PLA_free( senddispls );

      break;
    }
  case PLA_PROJ_ONTO_ROW:
    /* Case: matrix is marked projected onto rows. */

    PLA_Obj_global_align_col( Obj_from, &align_col_from );
    PLA_Obj_global_align_row( Obj_to,   &align_row_to );
    if ( align_col_from == align_row_to ){
      /* Case: matrix and multivector a perfectly aligned. */

      /* Get global and local width of the objects */
      PLA_Obj_global_length( Obj_to, &length_mv );
      PLA_Obj_local_length( Obj_to, &local_length_mv );

      /* Extract me, myrow, nprows, npcols, and col communicator */
      PLA_Obj_template( Obj_from, &templ );
      PLA_Temp_comm_col_rank( templ, &myrow );  /* myrow = row index in mesh */
      PLA_Temp_comm_all_rank( templ, &me );     /* me = index in all nodes   */
      PLA_Temp_comm_col( templ, &comm );        /* col communicator          */
      PLA_Temp_comm_col_size( templ, &nprows ); /* number of rows    in mesh */
      PLA_Temp_comm_row_size( templ, &npcols ); /* number of columns in mesh */

      /* ldims[ i ] will hold the number of cols of the target multivector
         that node ( i, mycol ) holds.                                       */

      ldims      = (int *) PLA_calloc( nprows, sizeof(int) ); 

      /* sendcounts[ i ] will hold the number of bytes to be sent to 
         node ( i, mycol ) ( = ldims[i] * size, where size will be 
         the width of the current panel being redistributed)                 */

      sendcounts = (int *) PLA_calloc( nprows, sizeof(int) );

      /* senddispls[ i ] will hold the displacement where the data to be
         sent to node ( i, mycol ) starts                                    */

      senddispls = (int *) PLA_calloc( nprows, sizeof(int) );

      /* Determine the index (relative to all nodes) of the first
         node in this node's col of nodes                                    */

      mycol_first_rank = me/nprows * nprows;

      /* compute ldims[i] = the local width of the target multivector on
         node ( i, mycol )                                                   */

      for ( i=0; i<nprows; i++ )
        PLA_Temp_vector_distant_length(templ,
                                  mycol_first_rank+i, 
                                  length_mv,
                                  align_row_to, 
                                  &ldims[i]);

      /* Take a views into Obj_from and Obj_to so that we can march
         through them slicing off row panels of Obj_from that all 
         exist within one row of nodes                                    */

      PLA_Obj_view_all( Obj_from, &A_cur );
      PLA_Obj_view_all( Obj_to,   &mv_cur );

      while( TRUE ){
	/* Determine length of row panel that exists all on one row of
           nodes                                                             */

	PLA_Obj_split_size( A_cur, PLA_SIDE_TOP, &size, &owner_row ); 

	if ( 0 == size ) break;

	/* Split off current panel to be redistributed                       */

	PLA_Obj_horz_split_2( A_cur,  size,    &A_1,  
                                                &A_cur );
	PLA_Obj_vert_split_2( mv_cur, size,    &mv_1, &mv_cur );

	/* Create a workspace into which to pack the current panel */
        PLA_Matrix_create_conf_to( A_1, &A_1_copy );

	/* Compute the sendcounts and displacements for the scatter 
           operation                                                         */
	senddispls[0] = 0;
	for ( i=0; i<nprows; i++ ) {
	  sendcounts[i] = ldims[i] * size;
	  if ( i>0 ) senddispls[i] = senddispls[i-1] + sendcounts[i-1];
	}

	/* If this node owns part of current panel, pack                     */

	if ( myrow == owner_row ){
	  PLA_Obj_local_buffer( A_1_copy, &sendbuf );

	  PLA_Pack_from_pmv_onto_row_to_mv( 
             A_1, sendbuf, sendcounts, senddispls, mv_1 );
	}

	/* Create a workspace into which to receive this node's 
           part of the multivector */

	PLA_Mvector_create_conf_to( mv_1, size, &mv_1_copy );
	PLA_Obj_local_buffer( mv_1_copy, &recvbuf );

	/* Scatter the data within cols.                                     */

	MPI_Scatterv( BF( sendbuf ), sendcounts, senddispls, datatype, 
                      BF( recvbuf ), local_length_mv * size, datatype, 
                      owner_row, comm );

	/* Copy into the multivector.  (Notice: we cannot receive 
           straight into mv_1 since the local leading dimension
           may not equal the number of rows being received.)                 */
        PLA_Unpack_from_pmv_onto_row_to_mv( mv_1, recvbuf );
/*	PLA_Local_copy( mv_1_copy, mv_1 ); */
      }

      /* Free temporary objects                                              */
      PLA_Obj_free( &A_cur );        
      PLA_Obj_free( &A_1 );
      PLA_Obj_free( &mv_cur );       
      PLA_Obj_free( &mv_1 );
      PLA_Obj_free( &A_1_copy );     
      PLA_Obj_free( &mv_1_copy );

      /* Free temporary arrays                                               */
      PLA_free( ldims );
      PLA_free( sendcounts );
      PLA_free( senddispls );

      break;
    }
  default:
    if ( report_case ) {
      printf("This case of PLA_Copy_from_matrix_to_mv not yet implemented\n");
      report_case = FALSE;
    }
    return PLA_Copy_old( Obj_from, Obj_to );
  }

}



/*----------------------------------------------------------------------*/

int PLA_Copy_from_matrix_to_pmv( PLA_Obj Obj_from, PLA_Obj Obj_to )

/************************************************************************
 
  Copy from matrix to projected multivector

*************************************************************************/

{
  if ( report_case ) {
    printf("PLA_Copy_from_matrix_to_pmv not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}

/*----------------------------------------------------------------------*/

int PLA_Copy_from_msc( PLA_Obj Obj_from, PLA_Obj Obj_to )

/************************************************************************
 
  Copy from multiscalar

*************************************************************************/

{
  int objtype_to, proj_onto, owner_col, owner_row;

  PLA_Obj_objtype( Obj_to, &objtype_to );
  
  switch( objtype_to ){
  case PLA_MSCALAR: 
    PLA_Copy_from_msc_to_msc( Obj_from, Obj_to );
    break;
/*  case PLA_VECTOR:   */
  case PLA_PMVECTOR: 
    PLA_Obj_project_onto( Obj_to, &proj_onto );
    switch( proj_onto ){
    case PLA_PROJ_ONTO_ROW: 
    case PLA_PROJ_ONTO_COL:
      PLA_Obj_owner_col( Obj_to, &owner_col );
      PLA_Obj_owner_row( Obj_to, &owner_row );
      if ( owner_col == PLA_ALL_COLS ||
	  owner_row == PLA_ALL_ROWS )
	PLA_Copy_from_msc_to_dpmv( Obj_from, Obj_to );
      else
	PLA_Copy_from_msc_to_pmv( Obj_from, Obj_to );
      break;
    }
    break;
  case PLA_MVECTOR: 
    PLA_Copy_from_msc_to_mv( Obj_from, Obj_to );
    break;
  case PLA_MATRIX:  
    PLA_Copy_from_msc_to_matrix( Obj_from, Obj_to );
    break;
  default:     printf("Obj_to in PLA_Copy has unknown objtype\n");
               exit( 0 );
  }

  return( PLA_SUCCESS );
}



/*----------------------------------------------------------------------*/

int PLA_Copy_from_msc_to_dpmv( PLA_Obj Obj_from, PLA_Obj Obj_to )

/************************************************************************
 
  Copy from multiscalar to duplicated projected multivector

*************************************************************************/

{
  if ( report_case ) {
    printf("PLA_Copy_from_msc_to_dpmv not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}


int PLA_Copy_from_matrix_to_dpmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  int
    proj_onto_from, proj_onto_to,
    align_from, align_to,
    local_length, local_width, local_ldim_to,
    myrow, mycol, owner,
    size, done = FALSE,
    typesize;
  PLA_Obj
    Obj_from_cur = NULL, Obj_from_1 = NULL,
    Obj_to_cur   = NULL, Obj_to_1   = NULL,
    Obj_temp = NULL;
  MPI_Datatype
    datatype;
  PLA_Template
    templ = NULL;
  MPI_Comm
    comm = MPI_COMM_NULL;
  void 
    *buffer_to, *temp_buffer;

  PLA_Obj_project_onto( Obj_from, &proj_onto_from );
  PLA_Obj_project_onto( Obj_to,   &proj_onto_to );

  if ( proj_onto_from == proj_onto_to ){
    if ( proj_onto_from == PLA_PROJ_ONTO_COL ){
      PLA_Obj_global_align_row( Obj_from, &align_from );
      PLA_Obj_global_align_row( Obj_to, &align_to );

      if ( align_from == align_to ){
	PLA_Obj_datatype( Obj_from, &datatype );
	MPI_Type_size( datatype, &typesize );

	PLA_Obj_template( Obj_from, &templ );
	PLA_Temp_comm_row( templ, &comm );
	PLA_Temp_comm_row_rank( templ, &mycol );

	PLA_Obj_local_length( Obj_from, &local_length );
	PLA_Obj_local_ldim  ( Obj_to, &local_ldim_to );
	PLA_Obj_view_all( Obj_from, &Obj_from_cur );
	PLA_Obj_view_all( Obj_to,   &Obj_to_cur );
	while ( TRUE ){
	  PLA_Obj_split_size_to_next_proc( Obj_from_cur, PLA_SIDE_LEFT, &size, &owner );
	  if ( 0 == size ) break;

	  PLA_Obj_vert_split_2( Obj_from_cur, size, 
                                          &Obj_from_1, &Obj_from_cur );
	  PLA_Obj_vert_split_2( Obj_to_cur, size, 
                                          &Obj_to_1,   &Obj_to_cur );
	  if ( mycol == owner ){
/*	    PLA_Local_copy( Obj_from_1, Obj_to_1 ); */
	    PLA_Obj_local_buffer( Obj_to_1, &buffer_to );
	    PLA_Obj_get_local_contents( Obj_from_1, PLA_NO_TRANS,
				       &local_length, &size, buffer_to, 
				       local_ldim_to, 1 );
	  }

	  PLA_Obj_local_buffer( Obj_to_1, &buffer_to );
	  if ( local_length == local_ldim_to ){
	    MPI_Bcast( BF( buffer_to ), size * local_length, datatype, 
		       owner, comm );
	  }
	  else{
	    temp_buffer = PLA_malloc( (size_t) typesize * size * local_length );
	    if ( mycol == owner ) 
	      PLA_Obj_get_local_contents( Obj_to_1, PLA_NO_TRANS,
					  &local_length, &size, temp_buffer,
					  local_length, 1 );
	    
	    MPI_Bcast( BF( temp_buffer ), size * local_length, datatype, 
		       owner, comm );
	    
	    if ( mycol != owner )
	      PLA_Obj_set_local_contents( PLA_NO_TRANS, 
					  local_length, size, temp_buffer,
					  local_length, 1, Obj_to_1 );

	    PLA_free( temp_buffer );
	  }
	}
	PLA_Obj_free( &Obj_from_1 );
	PLA_Obj_free( &Obj_from_cur );
	PLA_Obj_free( &Obj_to_1 );
	PLA_Obj_free( &Obj_to_cur );

	done = TRUE;
      }
    }
    else { /* proj_onto_from == PLA_PROJ_ONTO_ROW */ 
      PLA_Obj_global_align_col( Obj_from, &align_from );
      PLA_Obj_global_align_col( Obj_to, &align_to );

      if ( align_from == align_to ){
	PLA_Obj_datatype( Obj_from, &datatype );
	MPI_Type_size( datatype, &typesize );

	PLA_Obj_template( Obj_from, &templ );
	PLA_Temp_comm_col( templ, &comm );
	PLA_Temp_comm_col_rank( templ, &myrow );

	PLA_Obj_local_width( Obj_from, &local_width );
	PLA_Obj_local_ldim  ( Obj_to, &local_ldim_to );

	PLA_Obj_view_all( Obj_from, &Obj_from_cur );
	PLA_Obj_view_all( Obj_to,   &Obj_to_cur );

	while ( TRUE ){
	  PLA_Obj_split_size( Obj_from_cur, PLA_SIDE_TOP, &size, &owner );
	  if ( 0 == size ) break;

	  PLA_Obj_horz_split_2( Obj_from_cur, size, &Obj_from_1, 
                                                    &Obj_from_cur );

	  PLA_Obj_horz_split_2( Obj_to_cur, size,   &Obj_to_1, 
                                                    &Obj_to_cur );
	  PLA_Obj_local_buffer( Obj_to_1, &buffer_to );
	  if ( myrow == owner ){
/*	    PLA_Local_copy( Obj_from_1, Obj_to_1 ); */
	    PLA_Obj_get_local_contents( Obj_from_1, PLA_NO_TRANS,
				       &size, &local_width, buffer_to, 
				       local_ldim_to, 1 );
	  }

	  if ( size == local_ldim_to ){
	    MPI_Bcast( BF( buffer_to ), size * local_width, datatype, 
		       owner, comm );
	  }
	  else{
	    temp_buffer = PLA_malloc( (size_t) typesize * size * local_width );
	    if ( myrow == owner ) 
	      PLA_Obj_get_local_contents( Obj_to_1, PLA_NO_TRANS,
					  &size, &local_width, temp_buffer,
					  size, 1 );
	    
	    MPI_Bcast( BF( temp_buffer ), size * local_width, datatype, 
		       owner, comm );
	    
	    if ( myrow != owner )
	      PLA_Obj_set_local_contents( PLA_NO_TRANS, 
					  size, local_width, temp_buffer,
					  size, 1, Obj_to_1 );

	    PLA_free( temp_buffer );
	  }
	}
	PLA_Obj_free( &Obj_from_1 );
	PLA_Obj_free( &Obj_from_cur );
	PLA_Obj_free( &Obj_to_1 );
	PLA_Obj_free( &Obj_to_cur );

	done = TRUE;
      }
    }
  }
  else {
    if ( proj_onto_from == PLA_PROJ_ONTO_COL )
      PLA_Obj_global_width( Obj_from, &size );
    else
      PLA_Obj_global_length( Obj_from, &size );

    PLA_Mvector_create_conf_to( Obj_from, size, &Obj_temp );

    PLA_Copy_from_matrix_to_mv( Obj_from, Obj_temp );
    PLA_Copy_from_mv_to_dpmv( Obj_temp, Obj_to );
  
    PLA_Obj_free( &Obj_temp );
    done = TRUE; 
  }
  
  if ( !done ){
    if ( report_case ) {
      printf("PLA_Copy_from_matrix_to_dpmv not yet implemented\n");
      report_case = FALSE;
    }
    return PLA_Copy_old( Obj_from, Obj_to );
  }
  else
    return PLA_SUCCESS;
}



/*----------------------------------------------------------------------*/

int PLA_Copy_from_msc_to_matrix( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from multiscalar to matrix

*************************************************************************/
  int
    done = TRUE,
    owner_row, owner_col;

  PLA_Obj_owner_row( Obj_from, &owner_row );
  PLA_Obj_owner_col( Obj_from, &owner_col );

  if ( owner_row == PLA_ALL_ROWS && owner_col == PLA_ALL_COLS ){
    PLA_Copy_local_part_from_msc_to_matrix( Obj_from, Obj_to );
    
    done = TRUE;
  }
  else{
    PLA_Obj
      Obj_temp = NULL;

    PLA_Mscalar_create_conf_to( Obj_from, PLA_ALL_ROWS, PLA_ALL_COLS,
				 &Obj_temp );
    PLA_Copy_from_msc_to_msc( Obj_from, Obj_temp );
    PLA_Copy_from_msc_to_matrix( Obj_temp, Obj_to );
    
    PLA_Obj_free( &Obj_temp );
    
    done = TRUE;
  }
    

  if ( !done ){
    if ( report_case ) {
      printf("PLA_Copy_from_msc_to_matrix not yet implemented\n");
      report_case = FALSE;
    }
    return PLA_Copy_old( Obj_from, Obj_to );
  }
  
  return PLA_SUCCESS;
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_msc_to_msc( PLA_Obj Obj_from, PLA_Obj Obj_to )

/************************************************************************
 
  Copy from multiscalar to multiscalar

*************************************************************************/

{
  int
    done = FALSE,
    owner_row_from, owner_col_from,
    owner_row_to, owner_col_to,
    local_length, local_width, local_ldim,
    myrow, mycol, nprows, root;
  
  PLA_Template
    templ;

  void 
    *local_buffer;

  MPI_Datatype
    datatype;

  MPI_Comm
    comm = MPI_COMM_NULL;

  PLA_Obj_owner_row( Obj_from, &owner_row_from );
  PLA_Obj_owner_col( Obj_from, &owner_col_from );
  PLA_Obj_owner_row( Obj_to, &owner_row_to );
  PLA_Obj_owner_col( Obj_to, &owner_col_to );

  if ( owner_row_from == PLA_ALL_ROWS ) {
    if ( owner_col_from == PLA_ALL_COLS || owner_col_from == owner_col_to ) {
      PLA_Local_copy( Obj_from, Obj_to );
      done = TRUE;
    }
    else if ( owner_col_to == PLA_ALL_COLS ){
    }
    else{
    }
  }
  else if ( owner_row_from == owner_row_to ){
      if ( owner_col_from == PLA_ALL_COLS || owner_col_from == owner_col_to ){
	PLA_Local_copy( Obj_from, Obj_to );
	done = TRUE;
      }
      else if ( owner_col_to == PLA_ALL_COLS ){
      }
      else{
      }      
    }
  else{
    if ( owner_row_to == PLA_ALL_ROWS ){
      if ( owner_col_to == PLA_ALL_COLS ){
	PLA_Obj_datatype( Obj_from, &datatype );
	PLA_Obj_template( Obj_from, &templ );
	PLA_Temp_comm_col_rank( templ, &myrow );
	PLA_Temp_comm_col_size( templ, &nprows );
	PLA_Temp_comm_row_rank( templ, &mycol );
	PLA_Temp_comm_all( templ, &comm );
	if ( myrow == owner_row_from && mycol == owner_col_from )
	  PLA_Local_copy( Obj_from, Obj_to );
	PLA_Obj_local_length( Obj_to, &local_length );
	PLA_Obj_local_width( Obj_to, &local_width );
	PLA_Obj_local_ldim( Obj_to, &local_ldim );
	root = owner_col_from * nprows + owner_row_from;

	if ( local_length != local_ldim ){
	  MPI_Datatype
	    datatype;
	  int 
	    typesize;
	  void 
	    *buff_temp;

	  PLA_Obj_datatype( Obj_to, &datatype );
	  MPI_Type_size( datatype, &typesize );
	  buff_temp = PLA_malloc( (size_t) typesize * local_length * 
                                              local_width );
	  
	  if ( myrow == owner_row_from && mycol == owner_col_from ){
	    PLA_Obj_get_local_contents( Obj_to, PLA_NO_TRANS,
					&local_length, &local_width, 
					buff_temp, local_length, 1 );
	    MPI_Bcast( BF( buff_temp ), local_length * local_width,
		       datatype, root, comm );
	  }
	  else{
	    MPI_Bcast( BF( buff_temp ), local_length * local_width,
		       datatype, root, comm );

	    PLA_Obj_set_local_contents( PLA_NO_TRANS, 
					 local_length, local_width, 
                                         buff_temp, local_length, 1,
					 Obj_to );
	  }
	  PLA_free( buff_temp );

	  done = TRUE;
	}
	else {
	  PLA_Obj_local_buffer( Obj_to, &local_buffer );
	  MPI_Bcast( BF( local_buffer ), local_length*local_width,
		     datatype, root, comm );
	  done = TRUE;
	}
      }
    }
  }

  if ( !done ){
    if ( report_case ) {
      printf("PLA_Copy_from_msc_to_msc not yet implemented\n");
      report_case = FALSE;
    }

    return PLA_Copy_old( Obj_from, Obj_to );
  }
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_msc_to_mv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from multiscalar to multivector

*************************************************************************/


  int owner_row, owner_col;

  PLA_Obj_owner_row( Obj_from, &owner_row );
  PLA_Obj_owner_col( Obj_from, &owner_col );
  if ( owner_row == PLA_ALL_ROWS && owner_col == PLA_ALL_COLS ){
    int size, dummy;
    PLA_Obj   from_cur = NULL,        from_1 = NULL,
              to_cur   = NULL,        to_1   = NULL;    

    PLA_Obj_view_all( Obj_to,     &to_cur );
    PLA_Obj_view_all( Obj_from,   &from_cur );
  
    while ( TRUE ){
      PLA_Obj_split_size( to_cur, PLA_SIDE_TOP, &size, &dummy );
      if ( 0 == size ) break;
      PLA_Obj_horz_split_2( to_cur,   size,    &to_1,
                                               &to_cur );
      PLA_Obj_horz_split_2( from_cur, size,    &from_1,
                                               &from_cur );
      PLA_Obj_local_length( to_1, &size );
      if ( 0 != size ) PLA_Local_copy( from_1, to_1 );
    }

    PLA_Obj_free( &to_cur );      PLA_Obj_free( &to_1 );
    PLA_Obj_free( &from_cur );    PLA_Obj_free( &from_1 );
    
    return PLA_SUCCESS;
  }
  else if ( report_case ) {
    printf("PLA_Copy_from_msc_to_mv not yet implemented\n");
    report_case = FALSE;

    return PLA_Copy_old( Obj_from, Obj_to );
  }
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_msc_to_pmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from multiscalar to projected multivector

*************************************************************************/


  if ( report_case ) {
    printf("PLA_Copy_from_msc_to_ not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}

/*----------------------------------------------------------------------*/

int PLA_Copy_from_mv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from multivector

*************************************************************************/


  int objtype_to, proj_onto, owner_col, owner_row;

  PLA_Obj_objtype( Obj_to, &objtype_to );
  
  switch( objtype_to ){
  case PLA_MSCALAR: 
    PLA_Copy_from_mv_to_msc( Obj_from, Obj_to );
    break;
/*  case PLA_VECTOR:  */
  case PLA_PMVECTOR: 
    PLA_Obj_project_onto( Obj_to, &proj_onto );
    switch( proj_onto ){
    case PLA_PROJ_ONTO_ROW: 
    case PLA_PROJ_ONTO_COL:
      PLA_Obj_owner_col( Obj_to, &owner_col );
      PLA_Obj_owner_row( Obj_to, &owner_row );
      if ( owner_col == PLA_ALL_COLS ||
	  owner_row == PLA_ALL_ROWS )
	PLA_Copy_from_mv_to_dpmv( Obj_from, Obj_to );
      else
	PLA_Copy_from_mv_to_pmv( Obj_from, Obj_to );
    }
    break;
  case PLA_MVECTOR: 
    PLA_Copy_from_mv_to_mv( Obj_from, Obj_to );
    break;
  case PLA_MATRIX:  
    PLA_Copy_from_mv_to_matrix( Obj_from, Obj_to );
    break;
  default:
    printf("Obj_to in PLA_Copy has unknown objtype\n");
    exit( 0 );
  }

}





/*----------------------------------------------------------------------*/

int PLA_Copy_from_mv_to_dpmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from multivector to duplicated projected multivector

*************************************************************************/


  int 
    proj_onto_to, 
    align_row_from, align_row_to, align_col_to,
    local_length_mv, 
    width,
    sendcount,
    nprows, npcols, 
    *recvcounts, 
    *recvdispls;

  void 
     *sendbuf, 
     *recvbuf;
  PLA_Obj
    mv_copy = NULL,   
    dpmv_copy = NULL;

  PLA_Template 
    templ = NULL;

  MPI_Comm 
    comm = MPI_COMM_NULL;

  MPI_Datatype 
    datatype;

/*  PLA_Copy_old( Obj_from, Obj_to );
  return PLA_SUCCESS; */

  PLA_Obj_project_onto( Obj_to, &proj_onto_to );
  PLA_Obj_datatype( Obj_to, &datatype );

  PLA_Obj_template( Obj_from, &templ );

  switch( proj_onto_to ){
  case PLA_PROJ_ONTO_COL:
    PLA_Obj_global_align_row( Obj_from, &align_row_from );
    PLA_Obj_global_align_row( Obj_to,   &align_row_to );

    if ( align_row_from == align_row_to ){
      PLA_Obj_global_width( Obj_to, &width );

      PLA_Temp_comm_row_size( templ, &npcols );
      recvcounts = (int *) PLA_malloc( (size_t) npcols * sizeof(int) );
      recvdispls = (int *) PLA_malloc( (size_t) npcols * sizeof(int) );

      PLA_Compute_subbuffer_info_within_row( 
                          Obj_from, recvcounts, recvdispls );

      PLA_Mvector_create_conf_to( Obj_from, width, &mv_copy );
      PLA_Local_copy( Obj_from, mv_copy );
      PLA_Obj_local_buffer( mv_copy, &sendbuf );
      PLA_Obj_local_length( Obj_from, &local_length_mv );
      sendcount = local_length_mv * width;

      PLA_Pmvector_create_conf_to( Obj_to, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
				   width, &dpmv_copy );
      PLA_Obj_local_buffer( dpmv_copy, &recvbuf );

      PLA_Temp_comm_row( templ, &comm );
      
	  
      PLA_MPI_Allgatherv( BF(sendbuf) , sendcount, 
		   datatype, BF(recvbuf), recvcounts, 
		   recvdispls, datatype, comm ); 

      PLA_Unpack_from_mv_to_pmv_onto_col( 
             Obj_from, recvbuf, recvcounts, recvdispls, Obj_to );

      PLA_free( recvcounts );
      PLA_free( recvdispls );

      PLA_Obj_free( &mv_copy );
      PLA_Obj_free( &dpmv_copy );

      break;
    }
  case PLA_PROJ_ONTO_ROW:
    PLA_Obj_global_align_row( Obj_from, &align_row_from );
    PLA_Obj_global_align_col( Obj_to,   &align_col_to );

    if ( align_row_from == align_col_to ){
      PLA_Obj_global_width( Obj_from, &width );

      PLA_Temp_comm_col_size( templ, &nprows );
      recvcounts = (int *) PLA_malloc( (size_t) nprows * sizeof(int) );
      recvdispls = (int *) PLA_malloc( (size_t) nprows * sizeof(int) );

      PLA_Compute_subbuffer_info_within_col( 
                          Obj_from, recvcounts, recvdispls );

      PLA_Mvector_create_conf_to( Obj_from, width, &mv_copy );
      PLA_Obj_local_buffer( mv_copy, &sendbuf );
      PLA_Pack_from_mv_to_pmv_onto_row( Obj_from, sendbuf );
      PLA_Obj_local_length( Obj_from, &local_length_mv );
      sendcount = local_length_mv * width;

      PLA_Pmvector_create_conf_to( Obj_to, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
				   width, &dpmv_copy );
      PLA_Obj_local_buffer( dpmv_copy, &recvbuf );

      PLA_Temp_comm_col( templ, &comm );

      PLA_MPI_Allgatherv( BF(sendbuf), sendcount, 
		   datatype, BF(recvbuf), recvcounts, 
		   recvdispls, datatype, comm ); 


      PLA_Unpack_from_mv_to_pmv_onto_row( 
             Obj_from, recvbuf, recvcounts, recvdispls, Obj_to );

      PLA_free( recvcounts );
      PLA_free( recvdispls );

      PLA_Obj_free( &mv_copy );
      PLA_Obj_free( &dpmv_copy );
      break;
    }
  default:
    if ( report_case ) {
      printf("This case of PLA_Copy_from_mv_to_dpmv not yet implemented\n");
      report_case = FALSE;
    }
    return PLA_Copy_old( Obj_from, Obj_to );
  }
}



/*----------------------------------------------------------------------*/

int PLA_Copy_from_mv_to_matrix( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from multivector to matrix

  NOTICE: implementation suboptimal

*************************************************************************/
  int 
    done = FALSE,
    proj_onto,
    size;

  PLA_Obj
    Obj_temp = NULL;

  PLA_Obj_get_orientation( Obj_to, &proj_onto );
  PLA_Obj_global_width( Obj_from, &size );

  if ( proj_onto == PLA_PROJ_ONTO_COL )
    PLA_Pmvector_create_conf_to( Obj_to, PLA_PROJ_ONTO_COL, PLA_ALL_COLS,
				  size, &Obj_temp );
  else
    PLA_Pmvector_create_conf_to( Obj_to, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS,
				  size, &Obj_temp );

  PLA_Copy_from_mv_to_dpmv( Obj_from, Obj_temp );
  PLA_Copy_from_dpmv_to_matrix( Obj_temp, Obj_to );

  PLA_Obj_free( &Obj_temp );

  return PLA_SUCCESS;
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_mv_to_msc( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from multivector to multiscalar

*************************************************************************/
  int 
    done = FALSE,
    global_length,
    size_top,   
    owner_top;

  PLA_Obj_global_length( Obj_from, &global_length );
  PLA_Obj_split_size( Obj_from, PLA_SIDE_TOP, &size_top, &owner_top );

  if ( global_length == size_top ){
    /* Data exists on one node.  */
    PLA_Obj
      from_cur = NULL;

    PLA_Obj_view_all( Obj_from, &from_cur );
    PLA_Obj_objtype_cast( from_cur, PLA_MSCALAR );
    PLA_Copy_from_msc_to_msc( from_cur, Obj_to );

    PLA_Obj_free( &from_cur );

    done = TRUE;
  }
  else{
    PLA_Obj
      local_contr = NULL;

    PLA_Mscalar_create_conf_to( Obj_from, PLA_ALL_ROWS, PLA_ALL_COLS, 
			        &local_contr );

    PLA_Obj_set_to_zero( local_contr );
    PLA_Copy_local_part_from_mv_to_msc( Obj_from, local_contr );
    PLA_Reduce( local_contr, MPI_SUM, Obj_to );

    PLA_Obj_free( &local_contr );

    done = TRUE;
  }

  if ( !done ){
    if ( report_case ) {
      printf("PLA_Copy_from_mv_to_msc not yet implemented\n");
      report_case = FALSE;
    }
    return PLA_Copy_old( Obj_from, Obj_to );
  }

  return PLA_SUCCESS;
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_mv_to_mv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from multivector to multivector

*************************************************************************/


  if ( report_case ) {
    printf("PLA_Copy_from_mv_to_mv not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_mv_to_pmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from multivector to projected multivector

*************************************************************************/


  if ( report_case ) {
    printf("PLA_Copy_from_mv_to_pmv not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}

/*----------------------------------------------------------------------*/

int PLA_Copy_from_pmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from projected multivector

*************************************************************************/


  int objtype_to, proj_onto, owner_col, owner_row;

  PLA_Obj_objtype( Obj_to, &objtype_to );
  
  switch( objtype_to ){
  case PLA_MSCALAR: 
    PLA_Copy_from_pmv_to_msc( Obj_from, Obj_to );
    break;
/*  case PLA_VECTOR:  */
  case PLA_PMVECTOR: 
    PLA_Obj_project_onto( Obj_to, &proj_onto );
    switch( proj_onto ){
    case PLA_PROJ_ONTO_ROW: 
    case PLA_PROJ_ONTO_COL:
      PLA_Obj_owner_col( Obj_to, &owner_col );
      PLA_Obj_owner_row( Obj_to, &owner_row );
      if ( owner_col == PLA_ALL_COLS ||
	  owner_row == PLA_ALL_ROWS )
	PLA_Copy_from_pmv_to_dpmv( Obj_from, Obj_to );
      else
	PLA_Copy_from_pmv_to_pmv( Obj_from, Obj_to );
    }
    break;
  case PLA_MVECTOR: 
    PLA_Copy_from_pmv_to_mv( Obj_from, Obj_to );
    break;
  case PLA_MATRIX:  
    PLA_Copy_from_pmv_to_matrix( Obj_from, Obj_to );
    break;
  default:
    printf("Obj_to in PLA_Copy has unknown objtype\n");
    exit( 0 );
  }
}



/*----------------------------------------------------------------------*/

int PLA_Copy_from_pmv_to_dpmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from projected multivector to duplicated projected multivector

*************************************************************************/


  if ( report_case ) {
    printf("PLA_Copy_from_pmv_to_dpmv not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_pmv_to_matrix( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from projected multivector to matrix

*************************************************************************/


  if ( report_case ) {
    printf("PLA_Copy_from_pmv_to_matrix not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_pmv_to_msc( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from projected multivector to multiscalar

*************************************************************************/


  if ( report_case ) {
    printf("PLA_Copy_from_pmv_to_msc not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_pmv_to_mv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from projected multivector to multivector

*************************************************************************/


  if ( report_case ) {
    printf("PLA_Copy_from_pmv_to_mv not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}


/*----------------------------------------------------------------------*/

int PLA_Copy_from_pmv_to_pmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
/************************************************************************
 
  Copy from projected multivector to projected multivector

*************************************************************************/


  if ( report_case ) {
    printf("PLA_Copy_from_pmv_to_pmv not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}


int PLA_Copy_from_dpmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  int objtype_to, proj_onto, owner_col, owner_row;

  PLA_Obj_objtype( Obj_to, &objtype_to );
  
  switch( objtype_to ){
  case PLA_MSCALAR: 
    PLA_Copy_from_dpmv_to_msc( Obj_from, Obj_to );
    break;
/*  case PLA_VECTOR:  */
  case PLA_PMVECTOR: 
    PLA_Obj_project_onto( Obj_to, &proj_onto );
    switch( proj_onto ){
    case PLA_PROJ_ONTO_ROW: 
    case PLA_PROJ_ONTO_COL:
      PLA_Obj_owner_col( Obj_to, &owner_col );
      PLA_Obj_owner_row( Obj_to, &owner_row );
      if ( owner_col == PLA_ALL_COLS ||
	  owner_row == PLA_ALL_ROWS )
	PLA_Copy_from_dpmv_to_dpmv( Obj_from, Obj_to );
      else
	PLA_Copy_from_dpmv_to_pmv( Obj_from, Obj_to );
    }
    break;
  case PLA_MVECTOR: 
    PLA_Copy_from_dpmv_to_mv( Obj_from, Obj_to );
    break;
  case PLA_MATRIX:  
    PLA_Copy_from_dpmv_to_matrix( Obj_from, Obj_to );
    break;
  default:
    printf("Obj_to in PLA_Copy has unknown objtype\n");
    exit( 0 );
  }
}


int PLA_Copy_from_dpmv_to_msc( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  if ( report_case ) {
    printf("PLA_Copy_from_dpmv_to_msc not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}


int PLA_Copy_from_dpmv_to_matrix( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  int done = FALSE, proj_onto_from, proj_onto_to;

  PLA_Obj_project_onto( Obj_from, &proj_onto_from );
  PLA_Obj_project_onto( Obj_to, &proj_onto_to );
  if ( proj_onto_from == proj_onto_to ) {
    switch ( proj_onto_from ) {
    case PLA_PROJ_ONTO_COL:
      {
	int owner, mycol, size, align_row_from, align_row_to;
	PLA_Template templ = NULL;
	PLA_Obj A_dup_pmv_cur = NULL, A_mat_cur = NULL, A_dup_pmv_left = NULL, 
	A_mat_left = NULL;
 
	PLA_Obj_global_align_row( Obj_from, &align_row_from );
	PLA_Obj_global_align_row( Obj_to,   &align_row_to );
	if ( align_row_from != align_row_to ) break;
	else done = TRUE;

	PLA_Obj_template( Obj_from, &templ );
	PLA_Temp_comm_row_rank( templ, &mycol );

	PLA_Obj_view_all( Obj_from, &A_dup_pmv_cur );
	PLA_Obj_view_all( Obj_to, &A_mat_cur );
	while (TRUE ) {
	  PLA_Obj_split_size( A_mat_cur, PLA_SIDE_LEFT, &size, &owner );
	  if ( size == 0 ) break;
	  PLA_Obj_vert_split_2( A_mat_cur, size, &A_mat_left, &A_mat_cur );
	  PLA_Obj_vert_split_2( A_dup_pmv_cur, size, 
			       &A_dup_pmv_left, &A_dup_pmv_cur );
	  if ( owner == mycol ) PLA_Local_copy( A_dup_pmv_left, A_mat_left );
	}
	PLA_Obj_free( &A_dup_pmv_cur );
	PLA_Obj_free( &A_mat_cur );
	PLA_Obj_free( &A_dup_pmv_left );
	PLA_Obj_free( &A_mat_left );
	break;
      }
    case PLA_PROJ_ONTO_ROW:
      {
	int owner, myrow, size, align_col_from, align_col_to;
	PLA_Template templ = NULL;
	PLA_Obj A_dup_pmv_cur = NULL, A_mat_cur = NULL, A_dup_pmv_top = NULL, 
	A_mat_top = NULL;
 
	PLA_Obj_global_align_col( Obj_from, &align_col_from );
	PLA_Obj_global_align_col( Obj_to,   &align_col_to );
	if ( align_col_from != align_col_to ) break;
	else done = TRUE;

	PLA_Obj_template( Obj_from, &templ );
	PLA_Temp_comm_col_rank( templ, &myrow );

	PLA_Obj_view_all( Obj_from, &A_dup_pmv_cur );
	PLA_Obj_view_all( Obj_to, &A_mat_cur );
	while (TRUE ) {
	  PLA_Obj_split_size( A_mat_cur, PLA_SIDE_TOP, &size, &owner );
	  if ( size == 0 ) break;
	  PLA_Obj_horz_split_2( A_mat_cur, size, &A_mat_top, 
                                                 &A_mat_cur );
	  PLA_Obj_horz_split_2( A_dup_pmv_cur, size, &A_dup_pmv_top, 
                                                     &A_dup_pmv_cur );
	  if ( owner == myrow ) PLA_Local_copy( A_dup_pmv_top, A_mat_top );
	}
	PLA_Obj_free( &A_dup_pmv_cur );
	PLA_Obj_free( &A_mat_cur );
	PLA_Obj_free( &A_dup_pmv_top );
	PLA_Obj_free( &A_mat_top );
	break;
      }
    }
  }

  if ( !done ) {
    if ( report_case ) {
      printf("PLA_Copy_from_dpmv_to_matrix case not yet implemented\n");
      report_case = FALSE;
    }
    return PLA_Copy_old( Obj_from, Obj_to );
  }
  else return( PLA_SUCCESS );
}


int PLA_Copy_from_dpmv_to_mv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  int done = FALSE, proj_onto;

  PLA_Obj_project_onto( Obj_from, &proj_onto );

  switch( proj_onto ){
  case PLA_PROJ_ONTO_COL:
    {
      PLA_Obj to_cur = NULL,    from_cur = NULL,
              to_top = NULL,    from_top = NULL;
      int align_row_from, align_row_to,
          me, np, nb_distr, size, owner;
      PLA_Template templ = NULL;

      PLA_Obj_global_align_row( Obj_from, &align_row_from );
      PLA_Obj_global_align_row( Obj_to,   &align_row_to );
      if ( align_row_from == align_row_to ){
	PLA_Obj_template( Obj_from, &templ );
	PLA_Temp_comm_all_rank( templ, &me );
	PLA_Temp_comm_all_size( templ, &np );
	PLA_Temp_nb( templ, &nb_distr );

	PLA_Obj_view_all( Obj_to,   &to_cur );
	PLA_Obj_view_all( Obj_from, &from_cur );

	while( TRUE ){
	  PLA_Obj_split_size( to_cur, PLA_SIDE_TOP, &size, &owner );
	  if ( size == 0 ) break;
	  if ( me == owner ) {
	    PLA_Obj_horz_split_2( to_cur, size, &to_top,
				                &to_cur );
	    PLA_Obj_horz_split_2( from_cur, size, &from_top,
                                                  &from_cur );

	    PLA_Local_copy( from_top, to_top );

	    PLA_Obj_global_length( to_cur, &size );
	    if ( size < nb_distr * (np-1) ) break;

	    PLA_Obj_horz_split_2( to_cur, nb_distr*(np-1), PLA_DUMMY,
				                           &to_cur );
	    PLA_Obj_horz_split_2( from_cur, nb_distr*(np-1), PLA_DUMMY,
                                                             &from_cur );
	  }
	  else {
	    PLA_Obj_horz_split_2( to_cur, size, PLA_DUMMY,
				                &to_cur );
	    PLA_Obj_horz_split_2( from_cur, size, PLA_DUMMY,
                                                  &from_cur );
	  }
	}
	done = TRUE;

	PLA_Obj_free( &to_cur );             PLA_Obj_free( &from_cur );
	PLA_Obj_free( &to_top );             PLA_Obj_free( &from_top );
      }
      break;
    }
  case PLA_PROJ_ONTO_ROW:
    {
      PLA_Obj to_cur = NULL,    from_cur = NULL,
              to_top = NULL,    from_left = NULL;
      int align_col_from, align_row_to,
          me, np, nb_distr, size, owner, width_to, local_ldim_from;
      void *buffer_from;
      PLA_Template templ = NULL;

      PLA_Obj_global_align_col( Obj_from, &align_col_from );
      PLA_Obj_global_align_row( Obj_to,   &align_row_to );
      if ( align_col_from == align_row_to ){
	PLA_Obj_global_width( Obj_to, &width_to );
	PLA_Obj_local_ldim( Obj_from, &local_ldim_from );
	PLA_Obj_template( Obj_from, &templ );
	PLA_Temp_comm_all_rank( templ, &me );
	PLA_Temp_comm_all_size( templ, &np );
	PLA_Temp_nb( templ, &nb_distr );

	PLA_Obj_view_all( Obj_to,   &to_cur );
	PLA_Obj_view_all( Obj_from, &from_cur );

	while( TRUE ){
	  PLA_Obj_split_size( to_cur, PLA_SIDE_TOP, &size, &owner );
	  if ( size == 0 ) break;
	  if ( me == owner ) {
	    PLA_Obj_horz_split_2( to_cur, size, &to_top,
				                &to_cur );
	    PLA_Obj_vert_split_2( from_cur, size, &from_left, &from_cur );

/*	    PLA_Local_copy( from_top, to_top ); */
	    PLA_Obj_local_buffer( from_left, &buffer_from );
	    PLA_Obj_set_local_contents( PLA_TRANSPOSE, width_to, size,
				        buffer_from, local_ldim_from,
				        1, to_top );

	    PLA_Obj_global_length( to_cur, &size );
	    if ( size < nb_distr * (np-1) ) break;

	    PLA_Obj_horz_split_2( to_cur, nb_distr*(np-1), PLA_DUMMY,
				                           &to_cur );
	    PLA_Obj_vert_split_2( from_cur, nb_distr*(np-1), PLA_DUMMY,
                                                             &from_cur );
	  }
	  else {
	    PLA_Obj_horz_split_2( to_cur, size, PLA_DUMMY,
				                &to_cur );
	    PLA_Obj_vert_split_2( from_cur, size, PLA_DUMMY,
                                                  &from_cur );
	  }
	}
	done = TRUE;

	PLA_Obj_free( &to_cur );             PLA_Obj_free( &from_cur );
	PLA_Obj_free( &to_top );             PLA_Obj_free( &from_left );
      }
      break;
    }
  }

  if ( report_case ) {
    printf("PLA_Copy_from_dpmv_to_mv not yet implemented\n");
    report_case = FALSE;
  }
  if ( !done ) return PLA_Copy_old( Obj_from, Obj_to );
  else return( PLA_SUCCESS );
}


int PLA_Copy_from_dpmv_to_pmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  if ( report_case ) {
    printf("PLA_Copy_from_dpmv_to_pmv not yet implemented\n");
    report_case = FALSE;
  }
  return PLA_Copy_old( Obj_from, Obj_to );
}


int PLA_Copy_from_dpmv_to_dpmv( PLA_Obj Obj_from, PLA_Obj Obj_to )
{
  int done = FALSE, proj_onto_from, proj_onto_to;

  PLA_Obj_project_onto( Obj_from, &proj_onto_from );
  PLA_Obj_project_onto( Obj_to, &proj_onto_to );

  switch( proj_onto_from ){
  case PLA_PROJ_ONTO_COL:
    {
      if ( proj_onto_to == proj_onto_from ) {
	int align_row_from, align_row_to;

	PLA_Obj_global_align_row( Obj_from, &align_row_from );
	PLA_Obj_global_align_row( Obj_to,   &align_row_to );
	if ( align_row_from == align_row_to ){
	  PLA_Local_copy( Obj_from, Obj_to );
	  done = TRUE;
	  break;
	}
      }
      else {
	int align_row_from, align_col_to, width;

	PLA_Obj_global_align_row( Obj_from, &align_row_from );
	PLA_Obj_global_align_col( Obj_to,   &align_col_to );
	PLA_Obj_global_width( Obj_from, &width );
	if ( align_row_from == align_col_to ) {
	  PLA_Obj mv = NULL;

	  PLA_Mvector_create_conf_to( Obj_from, width, &mv );
	  PLA_Copy_from_dpmv_to_mv( Obj_from, mv );
	  PLA_Copy_from_mv_to_dpmv( mv, Obj_to );
	  done = TRUE;
	  
	  PLA_Obj_free( &mv );
	  break;
	}
      }
      break;
    }
  }
  if ( report_case ) {
    printf("PLA_Copy_from_dpmv_to_dpmv not yet implemented\n");
    report_case = FALSE;
  }
  if ( !done ) return PLA_Copy_old( Obj_from, Obj_to );
  else return( PLA_SUCCESS );
}



int PLA_Copy_old( PLA_Obj from, PLA_Obj to )
{
  int 
    objtype_from, objtype_to, owner_row, owner_col, owner;

  PLA_Obj_objtype( from, &objtype_from );
  PLA_Obj_objtype( to,   &objtype_to );

  switch( objtype_from ){
  case PLA_MATRIX:
    printf("Obj from: MATRIX\n");
    break;
  case PLA_PMVECTOR:
    printf("Obj from: PMVECTOR\n");
    break;
  case PLA_MVECTOR:
    printf("Obj from: MVECTOR\n");
    break;
  case PLA_MSCALAR:
    printf("Obj from: MSCALAR\n");

    PLA_Obj_owner_row ( from, &owner );
    if ( owner == PLA_ALL_ROWS ) 
      printf("owner_row = PLA_ALL_ROWS\n");
    else
      printf("owner_row = %d\n", owner);

    PLA_Obj_owner_col ( from, &owner );
    if ( owner == PLA_ALL_COLS ) 
      printf("owner_col = PLA_ALL_COLS\n");
    else
      printf("owner_col = %d\n", owner);

    break;
  default:
    printf("Unknown objtype_from\n");
  }

  switch( objtype_to ){
  case PLA_MATRIX:
    printf("Obj to: MATRIX\n");
    break;
  case PLA_PMVECTOR:
    printf("Obj to: PMVECTOR\n");
    break;
  case PLA_MVECTOR:
    printf("Obj to: MVECTOR\n");
    break;
  case PLA_MSCALAR:
    printf("Obj to: MSCALAR\n");

    PLA_Obj_owner_row ( to, &owner );
    if ( owner == PLA_ALL_ROWS ) 
      printf("owner_row = PLA_ALL_ROWS\n");
    else
      printf("owner_row = %d\n", owner);

    PLA_Obj_owner_col ( to, &owner );
    if ( owner == PLA_ALL_COLS ) 
      printf("owner_col = PLA_ALL_COLS\n");
    else
      printf("owner_col = %d\n", owner);

    break;
  default:
    printf("Unknown objtype_to\n");
  }

  PLA_Abort("PLA_Copy_old called: not implemented", __LINE__, __FILE__);
}

