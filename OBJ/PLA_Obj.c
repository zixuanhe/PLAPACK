/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/***************************************************************************/

int PLA_Vector_create   (
           MPI_Datatype   datatype,   int    global_length,
	   PLA_Template   templ,      int    global_align,
	   PLA_Obj       *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create distributed vector.


IN     datatype          datatype of object
IN     global_length     global length of vector
IN     templ             template for vector and matrix distribution
IN     global_align      alignment to template
OUT    new_obj           object describing created vector

----------------------------------------------------------------------------*/
{
  printf("Please use PLA_Mvector_create with width 1 to create vectors\n");
  
  exit( 0 );
}

/***************************************************************************/

int PLA_Mvector_create_without_buffer (
       MPI_Datatype   datatype,        int            global_length,
       int            global_width,    PLA_Template   templ,
       int            global_align,    PLA_Obj       *new_obj)
     
/*----------------------------------------------------------------------------

Purpose : Create distributed multi-vector, without buffer for data.

IN     datatype          datatype of object
IN     global_length     global length of multi-vector
IN     global_width      global width of multi-vector
IN     templ             template for vector and matrix distribution
IN     global_align      alignment to template
OUT    new_obj           object describing created multi-vector

----------------------------------------------------------------------------*/
{
  PLA_Obj_base 
    base_obj;

  PLA_Obj 
    view;

  int
    zero_or_one,
    value = PLA_SUCCESS;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Mvector_create_without_buffer_enter(
       datatype, global_length, global_width, templ, global_align, new_obj );

  /* if ( value == PLA_SUCCESS ) */{
    if ( *new_obj != NULL )
      PLA_Obj_free( new_obj );

    base_obj = pla_get_base_obj( );
    view     = pla_get_view( );
    
    base_obj->datatype         = datatype;

    base_obj->local_buffer     = NULL;

    base_obj->global_length    = global_length;
    base_obj->global_width     = global_width;

    PLA_Temp_zero_or_one( templ, &zero_or_one );
  
    base_obj->global_align_row = ( global_align == PLA_ALIGN_FIRST ? 
				   zero_or_one : global_align );
    base_obj->global_align_col = zero_or_one;
    base_obj->user_buffer      = PLA_UNDEFINED;
    base_obj->number_of_views  = 1;
    base_obj->local_buffer_size = 0;
    base_obj->ooc              = FALSE;
  
    base_obj->templ            = templ;

#if DEBUG==1
    view->pla_obj           = NULL;
#endif

    view->objtype           = PLA_MVECTOR;
    view->fixed             = FALSE; 
    view->global_length     = global_length;
    view->global_width      = global_width;
    view->global_align_row  = base_obj->global_align_row;
    view->global_align_col  = zero_or_one;
    view->owner_row         = PLA_ALL_ROWS;
    view->owner_col         = PLA_ALL_COLS;
    view->proj_onto         = PLA_PROJ_ONTO_COL;

    PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ALL,                   /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       view->global_align_row,        /* Alignment to template wrt direction */
       global_length,                 /* Length in direction */ 
       &view->local_length,           /* Length assigned to this node */
       &view->size_top,               /* Size of first block */
       &view->owner_top,	      /* Owner of first block */
       &view->size_bottom,            /* Size of last block */
       &view->owner_bottom );         /* Owner of last block */

    view->size_left         = global_width;
    view->owner_left        = 0;
    view->size_right        = global_width;
    view->owner_right       = 0;
    
    view->local_width       = global_width;
    base_obj->local_ldim    = view->local_length;

    view->local_stride      = 1;
  
    view->base_obj          = base_obj;

    view->mode   = PLA_MODE_CLOSED;

    *new_obj = view;
  }

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Mvector_create_without_buffer_exit(
       datatype, global_length, global_width, templ, global_align, new_obj );

  return value;
}

/***************************************************************************/

int PLA_Mvector_create   (
       MPI_Datatype   datatype,        int            global_length,
       int            global_width,    PLA_Template   templ,
       int            global_align,    PLA_Obj       *new_obj)

/*----------------------------------------------------------------------------

Purpose : Create distributed multi-vector.

IN     datatype          datatype of object
IN     global_length     global length of multi-vector
IN     global_width      global width of multi-vector
IN     templ             template for vector and matrix distribution
IN     global_align      alignment to template
OUT    new_obj           object describing created multi-vector

----------------------------------------------------------------------------*/
{
  int 
    value = PLA_SUCCESS;
  int local_length, local_width, type_size;
  void *buffer = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Mvector_create_enter(
       datatype, global_length, global_width, templ, global_align, new_obj );

  /* if ( value == PLA_SUCCESS ) */{
    PLA_Mvector_create_without_buffer( datatype, global_length, global_width,
				        templ, global_align, new_obj );
 
    PLA_OBJ_LOCAL_LENGTH( *new_obj, &local_length );
    PLA_OBJ_LOCAL_WIDTH ( *new_obj, &local_width );

    MPI_Type_size( datatype, &type_size );
    if ( local_length * local_width != 0 ){
      buffer = PLA_malloc( (size_t) type_size * local_length * local_width );
    
      if ( buffer == NULL ) 
	PLA_Abort( "PLA_Mvector_create: malloc failed\n", __LINE__, __FILE__);
    }

    PLA_Obj_attach_buffer( *new_obj, buffer, local_length, FALSE );

    PLA_Obj_set_to_zero( *new_obj );
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Mvector_create_exit(
       datatype, global_length, global_width, templ, global_align, new_obj );

  return value;
}

/***************************************************************************/

int PLA_Matrix_create_without_buffer   (
     MPI_Datatype   datatype,          int            global_length,
     int            global_width,      PLA_Template   templ,
     int            global_align_row,  int            global_align_col,
     PLA_Obj       *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create distributed matrix without buffer.

IN     datatype          datatype of object
IN     global_length     global length of matrix
IN     global_width      global width of matrix
IN     templ             template for vector and matrix distribution
IN     global_align_row  row alignment to template
IN     global_align_col  column alignment to template
OUT    new_obj           object describing created matrix

----------------------------------------------------------------------------*/
{
  PLA_Obj_base 
    base_obj;

  PLA_Obj 
    view;

  int
    zero_or_one,
    value = PLA_SUCCESS;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Matrix_create_without_buffer_enter(
       datatype, global_length, global_width, templ, 
       global_align_row, global_align_col, new_obj );

  /* if ( value == PLA_SUCCESS ) */{
    if ( *new_obj != NULL )
      PLA_Obj_free( new_obj );

    base_obj = pla_get_base_obj( );
    view     = pla_get_view( );
    
    base_obj->datatype         = datatype;

    base_obj->local_buffer     = NULL;

    base_obj->global_length    = global_length;
    base_obj->global_width     = global_width;

    PLA_Temp_zero_or_one( templ, &zero_or_one );
  
    base_obj->global_align_row = ( global_align_row == PLA_ALIGN_FIRST ? 
				   zero_or_one : global_align_row );
    base_obj->global_align_col = ( global_align_col == PLA_ALIGN_FIRST ? 
				   zero_or_one : global_align_col );
    base_obj->user_buffer      = PLA_UNDEFINED;
    base_obj->number_of_views  = 1;
    base_obj->local_buffer_size = 0;
    base_obj->ooc              = FALSE;
  
    base_obj->templ            = templ;

#if DEBUG==1
    view->pla_obj           = NULL;
#endif

    view->objtype           = PLA_MATRIX;
    view->fixed             = FALSE; 
    view->global_length     = global_length;
    view->global_width      = global_width;
    view->global_align_row  = base_obj->global_align_row;
    view->global_align_col  = base_obj->global_align_col;
    view->owner_row         = PLA_ALL_ROWS;
    view->owner_col         = PLA_ALL_COLS;
    view->proj_onto         = PLA_PROJ_ONTO_COL;

    PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_COL,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       view->global_align_row,        /* Alignment to template wrt direction */
       global_length,                 /* Length in direction */ 
       &view->local_length,           /* Length assigned to this node */
       &view->size_top,               /* Size of first block */
       &view->owner_top,	      /* Owner of first block */
       &view->size_bottom,            /* Size of last block */
       &view->owner_bottom );         /* Owner of last block */

    PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ROW,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       view->global_align_col,        /* Alignment to template wrt direction */
       global_width,                  /* Length in direction */ 
       &view->local_width,            /* Length assigned to this node */
       &view->size_left,              /* Size of first block */
       &view->owner_left,             /* Owner of first block */
       &view->size_right,             /* Size of last block */
       &view->owner_right );          /* Owner of last block */

  base_obj->local_ldim    = view->local_length;
  view->local_stride      = 1;
  
    view->base_obj          = base_obj;

    view->mode   = PLA_MODE_CLOSED;

    *new_obj = view;
  }

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Matrix_create_without_buffer_exit(
       datatype, global_length, global_width, templ, 
       global_align_row, global_align_col, new_obj );

  return value;
}

/***************************************************************************/

int PLA_Matrix_create   (
     MPI_Datatype   datatype,          int            global_length,
     int            global_width,      PLA_Template   templ,
     int            global_align_row,  int            global_align_col,
     PLA_Obj       *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create distributed matrix.

IN     datatype          datatype of object
IN     global_length     global length of matrix
IN     global_width      global width of matrix
IN     templ             template for vector and matrix distribution
IN     global_align_row  row alignment to template
IN     global_align_col  column alignment to template
OUT    new_obj           object describing created matrix

----------------------------------------------------------------------------*/
{
  int 
    value = PLA_SUCCESS;
  int local_length, local_width, type_size;
  void *buffer = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Matrix_create_enter(
       datatype, global_length, global_width, templ, 
       global_align_row, global_align_col, new_obj );

  /* if ( value == PLA_SUCCESS ) */{
    PLA_Matrix_create_without_buffer( datatype, global_length, global_width,
                      templ, global_align_row, global_align_col, new_obj );
  
    PLA_OBJ_LOCAL_LENGTH( *new_obj, &local_length );
    PLA_OBJ_LOCAL_WIDTH ( *new_obj, &local_width );

    MPI_Type_size( datatype, &type_size );
    if ( local_length * local_width != 0 ){
      buffer = PLA_malloc( (size_t) type_size * local_length * local_width );
    
      if ( buffer == NULL ) 
	PLA_Abort( "PLA_Matrix_create: malloc failed\n", __LINE__, __FILE__ );
    }

    PLA_Obj_attach_buffer( *new_obj, buffer, local_length, FALSE );

    PLA_Obj_set_to_zero( *new_obj );
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Matrix_create_exit(
       datatype, global_length, global_width, templ, 
       global_align_row, global_align_col, new_obj );

  return value;
}

/***************************************************************************/

int PLA_Mscalar_create_without_buffer  (
       MPI_Datatype   datatype,          int            owner_row,
       int            owner_col,         int            length,
       int            width,             PLA_Template   templ,
       PLA_Obj       *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create distributed multiscalar without buffer.


IN     datatype          datatype of object
IN     owner_row         index of row of nodes with contains owning node(s)
IN     owner_col         index of column of nodes with contains owning node(s)
IN     length            length of mscalar
IN     width             width of mscalar
IN     templ             template for vector and matrix distribution
OUT    new_obj           object describing created multiscalar

----------------------------------------------------------------------------*/
{
  PLA_Obj_base 
    base_obj;

  PLA_Obj 
    view;

  int
    zero_or_one,
    value = PLA_SUCCESS,
    myrow, mycol;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Mscalar_create_without_buffer_enter(
       datatype, owner_row, owner_col, length, width, templ, new_obj );

  /* if ( value == PLA_SUCCESS ) */{
    if ( *new_obj != NULL )
      PLA_Obj_free( new_obj );

    base_obj = pla_get_base_obj( );
    view     = pla_get_view( );
    
    base_obj->datatype         = datatype;

    base_obj->local_buffer     = NULL;

    base_obj->global_length    = length;
    base_obj->global_width     = width;

    base_obj->global_align_row = 0;
    base_obj->global_align_col = 0;

    base_obj->user_buffer      = PLA_UNDEFINED;
    base_obj->number_of_views  = 1;
    base_obj->local_buffer_size = 0;
    base_obj->ooc              = FALSE;
  
    base_obj->templ            = templ;

#if DEBUG==1
    view->pla_obj           = NULL;
#endif

    view->objtype           = PLA_MSCALAR;
    view->fixed             = FALSE; 
    view->global_length     = length;
    view->global_width      = width;
    view->global_align_row  = base_obj->global_align_row;
    view->global_align_col  = base_obj->global_align_col;
    view->owner_row         = owner_row;
    view->owner_col         = owner_col;
    view->proj_onto         = PLA_PROJ_ONTO_COL;

    PLA_Temp_comm_col_rank( templ, &myrow );
    PLA_Temp_comm_row_rank( templ, &mycol );

    /* Change this later! */
    if ( ( owner_row == PLA_ALL_ROWS || owner_row == myrow ) &&
         ( owner_col == PLA_ALL_COLS || owner_col == mycol ) ) {
      view->local_length = length; 
      view->local_width  = width;
      base_obj->local_ldim  = length;
    }
    else {
      view->local_length = 0;
      view->local_width  = 0;
      base_obj->local_ldim = 0;
    }

    view->size_top  = view->size_bottom = length;
    view->size_left = view->size_right  = width;

    view->owner_top  = view->owner_bottom = owner_row;
    view->owner_left = view->owner_right  = owner_col;

    view->local_stride      = 1;
  
    view->base_obj          = base_obj;

    view->mode   = PLA_MODE_CLOSED;

    *new_obj = view;
  }

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Mscalar_create_without_buffer_exit(
       datatype, owner_row, owner_col, length, width, templ, new_obj );

  return value;
}

/***************************************************************************/

int PLA_Mscalar_create  (
       MPI_Datatype   datatype,          int            owner_row,
       int            owner_col,         int            length,
       int            width,             PLA_Template   templ,
       PLA_Obj       *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create distributed mutltiscalar.


IN     datatype          datatype of object
IN     owner_row         index of row of nodes with contains owning node(s)
IN     owner_col         index of column of nodes with contains owning node(s)
IN     length            length of mscalar
IN     width             width of mscalar
IN     templ             template for vector and matrix distribution
OUT    new_obj           object describing created multiscalar

----------------------------------------------------------------------------*/
{
  int 
    value = PLA_SUCCESS;
  int local_length, local_width, type_size;
  void *buffer = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Mscalar_create_enter(
       datatype, owner_row, owner_col, length, width, templ, new_obj );

  /* if ( value == PLA_SUCCESS ) */{
    PLA_Mscalar_create_without_buffer( 
       datatype, owner_row, owner_col, length, width, templ, new_obj );
  
    PLA_OBJ_LOCAL_LENGTH( *new_obj, &local_length );
    PLA_OBJ_LOCAL_WIDTH ( *new_obj, &local_width );

    MPI_Type_size( datatype, &type_size );
    if ( local_length * local_width != 0 ){
      buffer = PLA_malloc( (size_t) type_size * local_length * local_width );
    
      if ( buffer == NULL ) 
	PLA_Abort( "PLA_Matrix_create: malloc failed\n", __LINE__, __FILE__ );
    }

    PLA_Obj_attach_buffer( *new_obj, buffer, local_length, FALSE );

    PLA_Obj_set_to_zero( *new_obj );
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Mscalar_create_exit(
       datatype, owner_row, owner_col, length, width, templ, new_obj );

  return value;
}

/***************************************************************************/

int PLA_Pvector_create  (
        MPI_Datatype   datatype,          int            project_onto,
        int            owner,             int            global_length,
        PLA_Template   templ,             int            global_align,
        PLA_Obj       *new_proj_vector)

/*----------------------------------------------------------------------------

Purpose : Create projected vector.


IN     datatype          datatype of object
IN     project_onto      direction onto which to project
IN     owner             index of row or column of nodes wich containts projected vector
IN     global_length     length of (projected) vector
IN     global_align      alignment to template
IN     templ             template for vector and matrix distribution
OUT    new_proj_vector   object describing created projected vector

----------------------------------------------------------------------------*/
{
  printf("Please use PLA_Pmvector_create with width 1 to create projected  vectors\n");
  
  exit( 0 );
}


/***************************************************************************/

int PLA_Pmvector_create_without_buffer  (
       MPI_Datatype   datatype,           int            project_onto,
       int            owner,              int            global_proj_length,
       int            global_proj_width,  PLA_Template   templ,
       int            global_align,       PLA_Obj       *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create projected multi-vector without buffer.

IN     datatype              datatype of object
IN     project_onto          direction onto which to project
IN     owner                 index of row or column of nodes wich containts projected vector
IN     global_proj_length    length of (projected) vector
IN     global_proj_width     width of (projected) vector
IN     global_align          alignment to template
IN     templ                 template for vector and matrix distribution
OUT    new_obj               object describing created projected multi-vector

----------------------------------------------------------------------------*/
{
  PLA_Obj_base 
    base_obj;

  PLA_Obj 
    view;

  int
    zero_or_one,
    value = PLA_SUCCESS,
    myrow, mycol;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Pmvector_create_without_buffer_enter(
       datatype, project_onto, owner, global_proj_length, global_proj_width, 
       templ, global_align, new_obj );

  /* if ( value == PLA_SUCCESS ) */{
    if ( *new_obj != NULL )
      PLA_Obj_free( new_obj );

    base_obj = pla_get_base_obj( );
    view     = pla_get_view( );
    
    base_obj->datatype         = datatype;

    base_obj->local_buffer     = NULL;

    base_obj->global_length    = global_proj_length;
    base_obj->global_width     = global_proj_width;

    PLA_Temp_zero_or_one( templ, &zero_or_one );
  
    if ( project_onto == PLA_PROJ_ONTO_COL ) {
      base_obj->global_align_row = ( global_align == PLA_ALIGN_FIRST ? 
                                     zero_or_one : global_align );
      base_obj->global_align_col = zero_or_one;
    }
    else {
      base_obj->global_align_col = ( global_align == PLA_ALIGN_FIRST ? 
                                     zero_or_one : global_align );
      base_obj->global_align_row = zero_or_one;
    }

    base_obj->user_buffer      = PLA_UNDEFINED;
    base_obj->number_of_views  = 1;
    base_obj->local_buffer_size = 0;
    base_obj->ooc              = FALSE;
  
    base_obj->templ            = templ;

#if DEBUG==1
    view->pla_obj           = NULL;
#endif

    view->objtype           = PLA_PMVECTOR;
    view->fixed             = FALSE; 
    view->global_length     = global_proj_length;
    view->global_width      = global_proj_width;
    view->global_align_row  = base_obj->global_align_row;
    view->global_align_col  = base_obj->global_align_col;
    view->proj_onto         = project_onto;
    if ( project_onto == PLA_PROJ_ONTO_COL ){
      view->owner_row       = PLA_ALL_ROWS;
      view->owner_col       = owner;
      PLA_Temp_compute_owners_and_sizes(
          PLA_COMM_COL,              /* Direction in which to compute */
          templ,                      /* Template for distribution */   
          view->global_align_row,     /* Alignment to template wrt direction */
          global_proj_length,         /* Length in direction */ 
          &view->local_length,        /* Length assigned to this node */
          &view->size_top,            /* Size of first block */
          &view->owner_top,	      /* Owner of first block */
          &view->size_bottom,         /* Size of last block */
          &view->owner_bottom );      /* Owner of last block */

      PLA_Temp_comm_row_rank( templ, &mycol );
      if ( mycol == owner || PLA_ALL_COLS == owner )
	view->local_width       = global_proj_width;
      else 
	view->local_width = 0;

      view->size_left         = global_proj_width;
      view->owner_left        = owner;
      view->size_right        = global_proj_width;
      view->owner_right       = owner;
    }
    else {
      view->owner_row       = owner;
      view->owner_col       = PLA_ALL_COLS;

      PLA_Temp_compute_owners_and_sizes(
          PLA_COMM_ROW,              /* Direction in which to compute */
          templ,                      /* Template for distribution */   
          view->global_align_col,     /* Alignment to template wrt direction */
          global_proj_width,          /* Length in direction */ 
          &view->local_width,         /* Length assigned to this node */
          &view->size_left,           /* Size of first block */
          &view->owner_left,	      /* Owner of first block */
          &view->size_right,          /* Size of last block */
          &view->owner_right );       /* Owner of last block */

      PLA_Temp_comm_col_rank( templ, &myrow );
      if ( myrow == owner || PLA_ALL_ROWS == owner )
	view->local_length   = global_proj_length;
      else 
	view->local_length = 0;

      view->size_top         = global_proj_length;
      view->owner_top        = owner;
      view->size_bottom      = global_proj_length;
      view->owner_bottom     = owner;
    }

    base_obj->local_ldim    = view->local_length;

    view->local_stride      = 1;
  
    view->base_obj          = base_obj;

    view->mode   = PLA_MODE_CLOSED;

    *new_obj = view;
  }

  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Pmvector_create_without_buffer_exit(
       datatype, project_onto, owner, global_proj_length, global_proj_width, 
       templ, global_align, new_obj );

  return value;
}

/***************************************************************************/

int PLA_Pmvector_create  (
       MPI_Datatype   datatype,           int            project_onto,
       int            owner,              int            global_proj_length,
       int            global_proj_width,  PLA_Template   templ,
       int            global_align,       PLA_Obj       *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create projected multi-vector.

IN     datatype              datatype of object
IN     project_onto          direction onto which to project
IN     owner                 index of row or column of nodes wich containts projected vector
IN     global_proj_length    length of (projected) vector
IN     global_proj_width     width of (projected) vector
IN     global_align          alignment to template
IN     templ                 template for vector and matrix distribution
OUT    new_obj               object describing created projected multi-vector

----------------------------------------------------------------------------*/
{
  int 
    value = PLA_SUCCESS;
  int local_length, local_width, type_size;
  void 
    *buffer = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Pmvector_create_enter(
       datatype, project_onto, owner, global_proj_length, global_proj_width, 
       templ, global_align, new_obj );

  /* if ( value == PLA_SUCCESS ) */{
    PLA_Pmvector_create_without_buffer( 
       datatype, project_onto, owner, global_proj_length, global_proj_width, 
       templ, global_align, new_obj );
  
    PLA_OBJ_LOCAL_LENGTH( *new_obj, &local_length );
    PLA_OBJ_LOCAL_WIDTH ( *new_obj, &local_width );

    MPI_Type_size( datatype, &type_size );
    if ( local_length * local_width != 0 ){
      buffer = PLA_malloc( (size_t) type_size * local_length * local_width );
    
      if ( buffer == NULL ) 
	PLA_Abort( "PLA_Matrix_create: malloc failed\n", __LINE__, __FILE__ );
    }

    PLA_Obj_attach_buffer( *new_obj, buffer, local_length, FALSE );

    PLA_Obj_set_to_zero( *new_obj );
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Pmvector_create_exit(
       datatype, project_onto, owner, global_proj_length, global_proj_width, 
       templ, global_align, new_obj );

  return value;
}

/***************************************************************************/

int PLA_Obj_free   (PLA_Obj   *obj)

/*----------------------------------------------------------------------------

Purpose : Linear algebra object destructor.


IN/OUT     object      object to be freed
-----------------------------------------------------------------------------*/
{
  int 
    value = PLA_SUCCESS;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_free_enter( obj );

  pla_free_view( *obj );

  *obj = NULL;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_free_exit( obj );

  return value;
}

/***************************************************************************/

int PLA_Obj_objtype  ( PLA_Obj     obj,    int *objtype)

/*----------------------------------------------------------------------------

Purpose : Extract object type from object.


IN         obj        object to be queried
OUT        objtype    object type of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_objtype_enter( obj, objtype );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_OBJTYPE( obj, objtype );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_objtype_exit( obj, objtype );

  return value;
}

/***************************************************************************/

int PLA_Obj_datatype  (PLA_Obj     obj,    MPI_Datatype *datatype)

/*----------------------------------------------------------------------------

Purpose : Extract datatype from object.


IN         obj         object to be queried
OUT        datatype    datatype of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_datatype_enter( obj, datatype );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_DATATYPE( obj, datatype );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_datatype_exit( obj, datatype );

  return value;
}

/***************************************************************************/

int PLA_Obj_template   (PLA_Obj   obj,     PLA_Template    *templ )

/*----------------------------------------------------------------------------

Purpose : Extract template from object.


IN         obj               object to be queried
OUT        templ             template

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_template_enter( obj, templ );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_TEMPLATE( obj, templ );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_template_exit( obj, templ );

  return value;
}

/***************************************************************************/

int PLA_Obj_global_info   (PLA_Obj   obj,             int    *global_length,
                           int       *global_width,   int    *project_onto,
                           int       *owner_row,      int    *owner_col,
                           int       *global_align_row,
                           int       *global_align_col)

/*----------------------------------------------------------------------------

Purpose :  Extract global information from linear algebra object.

Note that the individual fields can
be queried (separately) through calls to PLA_Obj_global(_field-name).


IN         obj                 object to be queried
OUT        global_length       global row dimension of object
OUT        global_width        global column dimension of object
OUT        project_onto        direction of projection
OUT        owner_row           row index of owning node(s)
OUT        owner_col           column index of owning node(s)
OUT        global_align(_row)  (row) alignment to template
OUT        global_align_col    column alignment to template

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_info_enter( obj, global_length, global_width,
                     project_onto, owner_row, owner_col, global_align_row,
                     global_align_col );

  /* if ( value == PLA_SUCCESS ) */{
    PLA_OBJ_GLOBAL_LENGTH( obj, global_length );
    PLA_OBJ_GLOBAL_WIDTH ( obj, global_width );
    PLA_OBJ_PROJECT_ONTO ( obj, project_onto );
    PLA_OBJ_OWNER_ROW    ( obj, owner_row );
    PLA_OBJ_OWNER_COL    ( obj, owner_col );
    PLA_OBJ_GLOBAL_ALIGN_ROW ( obj, global_align_row );
    PLA_OBJ_GLOBAL_ALIGN_COL ( obj, global_align_col );
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_info_exit( obj, global_length, global_width,
                     project_onto, owner_row, owner_col, global_align_row,
                     global_align_col );

  return value;
}

/***************************************************************************/

int PLA_Obj_global_length (PLA_Obj   obj,             int    *global_length )

/*----------------------------------------------------------------------------

Purpose :  Extract global length from linear algebra object.

IN         obj                 object to be queried
OUT        global_length       global row dimension of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_length_enter( obj, global_length );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_GLOBAL_LENGTH( obj, global_length );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_length_exit( obj, global_length );

  return value;
}

/***************************************************************************/

int PLA_Obj_global_width (PLA_Obj   obj,             int    *global_width )

/*----------------------------------------------------------------------------

Purpose :  Extract global width from linear algebra object.

IN         obj                 object to be queried
OUT        global_width        global column dimension of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_width_enter( obj, global_width );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_GLOBAL_WIDTH( obj, global_width );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_width_exit( obj, global_width );

  return value;
}

/***************************************************************************/

int PLA_Obj_global_size (PLA_Obj   obj,             int    *global_size )

/*----------------------------------------------------------------------------

Purpose :  Extract global size from linear algebra object.

IN         obj                 object to be queried
OUT        global_size         global row dimension of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  PLA_Abort( "PLA_Obj_global_size not supported.  Use global length or width.", __LINE__, __FILE__ );
}

/***************************************************************************/

int PLA_Obj_global_numvecs (PLA_Obj   obj,             int    *num_vecs )
{
  int
    value = PLA_SUCCESS;
  
  PLA_Abort( "PLA_Obj_global_numvecs not supported.  Use global length or width.", __LINE__, __FILE__ );
}

/***************************************************************************/

int PLA_Obj_project_onto  (PLA_Obj   obj,         int    *project_onto )

/*----------------------------------------------------------------------------

Purpose :  Extract projection information from linear algebra object.


IN         obj                 object to be queried
OUT        project_onto        direction of projection

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_project_onto_enter( obj, project_onto );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_PROJECT_ONTO( obj, project_onto );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_project_onto_exit( obj, project_onto );

  return value;
}

/***************************************************************************/

int PLA_Obj_owner_row     (PLA_Obj   obj,           int       *owner_row )

/*----------------------------------------------------------------------------

Purpose :  Extract mesh row owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_row           row index of owning node(s)
OUT        owner_col           column index of owning node(s)

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    objtype,
    size,
    global_length;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_owner_row_enter( obj, owner_row );

  /* if ( value == PLA_SUCCESS ) */ 
  PLA_Obj_objtype( obj, &objtype );

  if ( objtype == PLA_MATRIX ){
    PLA_Obj_split_size( obj, PLA_SIDE_TOP, &size, owner_row );
    PLA_Obj_global_length( obj, &global_length );
    if ( size != global_length )
      PLA_OBJ_OWNER_ROW( obj, owner_row );
  }
  else
    PLA_OBJ_OWNER_ROW( obj, owner_row );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_owner_row_exit( obj, owner_row );

  return value;
}

/***************************************************************************/

int PLA_Obj_owner_col     (PLA_Obj   obj,           int       *owner_col )

/*----------------------------------------------------------------------------

Purpose :  Extract mesh column owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_col           column index of owning node(s)

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_owner_col_enter( obj, owner_col );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_OWNER_COL( obj, owner_col );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_owner_col_exit( obj, owner_col );

  return value;
}

/***************************************************************************/

int PLA_Obj_global_align  (PLA_Obj   obj,      int       *global_align )

/*----------------------------------------------------------------------------

WARNING: check this out!!

Purpose :  Extract global alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align        Alignment to template

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_align_enter( obj, global_align );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_GLOBAL_ALIGN( obj, global_align );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_align_exit( obj, global_align );

  return value;
}

/***************************************************************************/

int PLA_Obj_global_align_row
                     (PLA_Obj   obj,      int       *global_align_row )

/*-------------------------------------------------------------------------

Purpose :  Extract global row alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_row  row alignment to template

-------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_align_row_enter( obj, global_align_row );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_GLOBAL_ALIGN_ROW( obj, global_align_row );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_align_row_exit( obj, global_align_row );

  return value;
}

/***************************************************************************/

int PLA_Obj_global_align_col
          (PLA_Obj   obj,      int       *global_align_col )

/*-------------------------------------------------------------------------

Purpose :  Extract global column alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_col    column alignment to template

-------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_align_col_enter( obj, global_align_col );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_GLOBAL_ALIGN_COL( obj, global_align_col );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_global_align_col_exit( obj, global_align_col );

  return value;
}

/***************************************************************************/

int PLA_Obj_local_info   (PLA_Obj   obj,             int    *local_length,
                          int       *local_width,    void   **local_buffer,
                          int       *local_stride,   int    *local_ldim)

/*----------------------------------------------------------------------------

Purpose :  Extract local information from linear algebra object.

Note that the individual fields can
be queried (separately) through calls to PLA_Obj_local(_field-name).



IN         obj                 object to be queried
OUT        local_length        local row dimension of object
OUT        local_width         local column dimension of object
OUT        local_buffer        address of local data
OUT        local_stride        stride between entries in a column or vector
OUT        local_ldim          leading dimension of array holding local data

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_info_enter( obj, local_length, local_width,
			       local_buffer, local_stride, local_ldim );

  /* if ( value == PLA_SUCCESS ) */ {
    PLA_OBJ_LOCAL_LENGTH( obj, local_length );
    PLA_OBJ_LOCAL_WIDTH ( obj, local_length );
    PLA_OBJ_LOCAL_BUFFER( obj, local_buffer );
    PLA_OBJ_LOCAL_STRIDE( obj, local_stride );
    PLA_OBJ_LOCAL_LDIM( obj, local_ldim );
  }
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_info_exit( obj, local_length, local_width,
			       local_buffer, local_stride, local_ldim );

  return value;
}

/***************************************************************************/

int PLA_Obj_local_length   (
	    PLA_Obj   obj,       int    *local_length )

/*----------------------------------------------------------------------------

Purpose :  Extract local length from linear algebra object.

IN         obj                 object to be queried
OUT        local_length        local row dimension of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_length_enter( obj, local_length );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_LOCAL_LENGTH( obj, local_length );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_length_exit( obj, local_length );

  return value;
}

/***************************************************************************/

int PLA_Obj_local_width   (PLA_Obj   obj,    int       *local_width )

/*----------------------------------------------------------------------------

Purpose :  Extract local width from linear algebra object.

IN         obj                 object to be queried
OUT        local_width         local column dimension of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_width_enter( obj, local_width );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_LOCAL_WIDTH( obj, local_width );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_width_exit( obj, local_width );

  return value;
}

/***************************************************************************/

int PLA_Obj_local_buffer  (PLA_Obj   obj,    void   **local_buffer )

/*----------------------------------------------------------------------------

Purpose :  Extract local buffer from linear algebra object.

IN         obj                 object to be queried
OUT        local_buffer        address of local data

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_buffer_enter( obj, local_buffer );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_LOCAL_BUFFER( obj, local_buffer );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_buffer_exit( obj, local_buffer );

  return value;
}

/***************************************************************************/

int PLA_Obj_local_ldim   (PLA_Obj   obj,             int    *local_ldim )

/*----------------------------------------------------------------------------

Purpose :  Extract local leading dimension from linear algebra object.

IN         obj                 object to be queried
OUT        local_ldim          leading dimension of array holding local data

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_ldim_enter( obj, local_ldim );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_LOCAL_LDIM( obj, local_ldim );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_ldim_exit( obj, local_ldim );

  return value;
  *local_ldim = obj->base_obj->local_ldim;

  return PLA_SUCCESS;
}

/***************************************************************************/

int PLA_Obj_local_stride   (PLA_Obj   obj,        int    *local_stride )

/*----------------------------------------------------------------------------

Purpose :  Extract local stride from linear algebra object.

IN         obj                 object to be queried
OUT        local_stride        stride between entries in a column or vector

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_stride_enter( obj, local_stride );

  /* if ( value == PLA_SUCCESS ) */ 
    PLA_OBJ_LOCAL_STRIDE( obj, local_stride );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_local_stride_exit( obj, local_stride );

  return value;
}

/***************************************************************************/

int PLA_Obj_get_local_contents   (
              PLA_Obj       obj,              int      trans,
              int            *rows_in_buf,     int      *cols_in_buf,
              void           *buf,             int      leading_dim_buf,
              int            stride_buf)

/*----------------------------------------------------------------------------

Purpose : Extract the local data from the given object.

IN     obj               global object
IN     trans             indicates whether to transpose data
OUT    rows_in_buf       row dimension of extracted data
OUT    cols_in_buf       column dimension of extracted data
OUT    buf               address where data is to be put
IN     leading_dim       leading dimension of buffer where data is put
IN     stride_buf        stride of buffer where data is put

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    m, n, i, j, k, typesize, ldim, stride;
    char *buf_local, *buf_obj, *tempp_local, *tempp_obj;
  MPI_Datatype 
    datatype;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_get_local_contents_enter( 
	        obj, trans, rows_in_buf, cols_in_buf, buf, 
                leading_dim_buf, stride_buf );

  /* if ( value == PLA_SUCCESS ) */ {
    PLA_Obj_local_length( obj, &m );
    PLA_Obj_local_width ( obj, &n );
    *rows_in_buf = m;
    *cols_in_buf = n;
    PLA_Obj_local_ldim  ( obj, &ldim );    
    PLA_Obj_datatype    ( obj, &datatype );
    MPI_Type_size( datatype, &typesize );
      
    buf_local = (char *) buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );
	
    if ( trans == PLA_NO_TRANS ){
      if ( leading_dim_buf == m && ldim == m )
	memcpy( buf_local, buf_obj, m*n*typesize ); 
      else{
	int
	  done = FALSE;

	if ( m == 1 ){
	  if ( datatype == MPI_DOUBLE ){
	    PLA_dcopy( &n, (double *) buf_obj, &ldim, (double *) buf, &leading_dim_buf );
	    done = TRUE;
	  }
	  else if ( datatype == MPI_FLOAT ){
	    PLA_scopy( &n, (float *) buf_obj, &ldim, (float *) buf, &leading_dim_buf );
	    done = TRUE;
	  }
	  else if ( datatype == MPI_DOUBLE_COMPLEX ){
	    PLA_zcopy( &n, (PLA_DOUBLE_COMPLEX *) buf_obj, &ldim, 
		       (PLA_DOUBLE_COMPLEX *) buf, &leading_dim_buf );
	    done = TRUE;
	  }
	  else if ( datatype == MPI_COMPLEX ){
	    PLA_ccopy( &n, ( PLA_COMPLEX * ) buf_obj, &ldim, 
		       ( PLA_COMPLEX * ) buf, &leading_dim_buf );
	    done = TRUE;
	  }
	}
	if ( !done ){
	  for ( j=0; j<n; j++ ){
	    tempp_local = buf_local + j*leading_dim_buf*typesize;
	    tempp_obj   = buf_obj   + j*ldim*typesize;
	    memcpy( tempp_local, tempp_obj, m*typesize ); 
	    /*	for ( i=0; i<m*typesize; i++ )
	     *tempp_local++ = *tempp_obj++;  */
	  }
	}
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = buf_local + j*leading_dim_buf*typesize;
	tempp_obj   = buf_obj + j*typesize;
	for ( i=0; i<n; i++ ) {
	  memcpy( tempp_local, tempp_obj, typesize ); 
	  tempp_obj += ldim*typesize; 
	  tempp_local += typesize; 
/*	  for ( k=0; k<typesize; k++ )
	    *tempp_local++ = *tempp_obj++;
	  tempp_obj += (ldim-1)*typesize;  */
	}
      }
    }
  }
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_get_local_contents_exit(
	        obj, trans, rows_in_buf, cols_in_buf, buf, 
                leading_dim_buf, stride_buf );

  return value;
}

/***************************************************************************/

int PLA_Obj_set_local_contents   (
              int      trans,                 int      rows_in_buf,
              int      cols_in_buf,           void     *buf,
              int      leading_dim_buf,       int      stride_buf,
              PLA_Obj obj )

/*----------------------------------------------------------------------------

Purpose : Set the local data from the given object.


IN        trans             indicates whether to transpose data
IN        rows_in_buf       row dimension of data buffer
IN        cols_in_buf       column dimension of data buffer
IN        buf               address of data buffer
IN        leading_dim       leading dimension of data buffer
IN        stride_buf        stride of data buffer where data is put
IN/OUT    obj               global object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    m, n, i, j, k, typesize, ldim, stride;
  void
    *buf_local, *buf_obj;
  char
    *tempp_local, *tempp_obj;
  MPI_Datatype 
    datatype;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_set_local_contents_enter( 
	        trans, rows_in_buf, cols_in_buf, buf, 
                leading_dim_buf, stride_buf, obj );

  /* if ( value == PLA_SUCCESS ) */ {
    PLA_Obj_local_length( obj, &m );
    PLA_Obj_local_width ( obj, &n );
    PLA_Obj_local_ldim  ( obj, &ldim );    
    PLA_Obj_datatype    ( obj, &datatype );
    MPI_Type_size( datatype, &typesize );
      
    buf_local = buf;
    PLA_Obj_local_buffer( obj, (void **) &buf_obj );

    if ( trans == PLA_NO_TRANS ){
      if ( leading_dim_buf == m && ldim == m )
	memcpy( (char *) buf_obj, (char *) buf_local, m*n*typesize ); 
      else{
	int
	  done = FALSE;

	if ( m == 1 ){
	  if ( datatype == MPI_DOUBLE ){
	    PLA_dcopy( &n, (double *) buf_local, &leading_dim_buf, 
		        (double *) buf_obj, &ldim );
	    done = TRUE;
	  }
	  else if ( datatype == MPI_FLOAT ){
	    PLA_scopy( &n, (float *) buf_local, &leading_dim_buf, 
		       (float *) buf_obj, &ldim );
	    done = TRUE;
	  }
	  else if ( datatype == MPI_DOUBLE_COMPLEX ){
	    PLA_zcopy( &n, (PLA_DOUBLE_COMPLEX *) buf_local, &leading_dim_buf, 
		       (PLA_DOUBLE_COMPLEX *) buf_obj, &ldim );
	    done = TRUE;
	  }
	  else if ( datatype == MPI_COMPLEX ){
	    PLA_ccopy( &n, (PLA_COMPLEX *) buf_local, &leading_dim_buf, 
		       (PLA_COMPLEX *) buf_obj, &ldim );
	    done = TRUE;
	  }
	}
	if ( !done ){
	  for ( j=0; j<n; j++ ){
	    tempp_local = (char *) buf_local + j*leading_dim_buf*typesize;
	    tempp_obj   = (char *) buf_obj   + j*ldim*typesize;
	    memcpy( tempp_obj, tempp_local, m*typesize ); 
	  }
	}
      }
    }
    else {
      for ( j=0; j<m; j++ ){
	tempp_local = (char *) buf_local + j*leading_dim_buf*typesize;
	tempp_obj   = (char *) buf_obj + j*typesize;
	for ( i=0; i<n; i++ ) {
	  memcpy( tempp_obj, tempp_local, typesize ); 
	  tempp_obj += ldim*typesize;
	  tempp_local += typesize;
/*	  for ( k=0; k<typesize; k++ )
	    *tempp_obj++ = *tempp_local++;
	  tempp_obj += (ldim-1)*typesize; */
	}
      }
    }
  }
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_set_local_contents_exit(
	        trans, rows_in_buf, cols_in_buf, buf, 
                leading_dim_buf, stride_buf, obj );

  return value;
}

/***************************************************************************/

int PLA_Obj_set   (PLA_Obj      obj,    MPI_Datatype   datatype,
                   void         *value )

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to value given
by value. The value is cast to the given MPI_Datatype via the datatype
parameter.


IN/OUT          obj          linear algebra object
IN              datatype     datatype of value
IN             *value        value to be used for initialization

----------------------------------------------------------------------------*/
{
  int
    return_value = PLA_SUCCESS,
    m, n, i, j, ldim;
  MPI_Datatype 
    datatype_obj;

  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_set_enter( obj, datatype, value );

  /* if ( return_value == PLA_SUCCESS ) */ {
    PLA_Obj_local_length( obj, &m );
    PLA_Obj_local_width ( obj, &n );
    PLA_Obj_local_ldim  ( obj, &ldim );    
    PLA_Obj_datatype    ( obj, &datatype_obj );

    if ( datatype_obj == MPI_DOUBLE ) {
      double 
	*buf_obj, *tempp, d_value;

      PLA_Obj_local_buffer( obj, (void **) &buf_obj );

      if ( datatype == MPI_DOUBLE )
	d_value = *( ( double *) value );
      else if ( datatype == MPI_FLOAT )
	d_value = ( double ) *( ( float *) value );
      else
	PLA_Abort( "Datatype not yet implemented", __LINE__, __FILE__ );

      for ( j=0; j<n; j++ ){
	tempp = buf_obj   + j*ldim;
	for ( i=0; i<m; i++ )
	  *tempp++ = d_value;
      }
    }
    else if ( datatype_obj == MPI_FLOAT ) {
      float
	*buf_obj, *tempp, f_value;

      PLA_Obj_local_buffer( obj, (void **) &buf_obj );

      if ( datatype == MPI_DOUBLE) 
	f_value = ( float ) *( ( double *) value );
      else if ( datatype == MPI_FLOAT )
	f_value = *( ( float *) value );
      else
	PLA_Abort( "Datatype not yet implemented", __LINE__, __FILE__ );
      
      for ( j=0; j<n; j++ ){
	tempp   = buf_obj   + j*ldim;
	for ( i=0; i<m; i++ )
	  *tempp++ = f_value;
      }
    }
    else if ( datatype_obj == MPI_DOUBLE_COMPLEX ) {
      double 
	*buf_obj, *tempp, d_value1, d_value2;

      PLA_Obj_local_buffer( obj, (void **) &buf_obj );
      
      if ( datatype == MPI_DOUBLE ){
	d_value1 = *( ( double *) value );
	d_value2 = 0.0;
      }
      else if ( datatype == MPI_FLOAT ){
	d_value1 = ( double ) *( ( float *) value );
	d_value2 = 0.0;
      }
      else 
	PLA_Abort( "Datatype not yet implemented", __LINE__, __FILE__ );

      for ( j=0; j<n; j++ ){
	tempp = buf_obj   + j*ldim*2;
	for ( i=0; i<m; i++ ){
	  *tempp++ = d_value1;
	  *tempp++ = d_value2;
	}
      }
    }
    else if ( datatype_obj == MPI_COMPLEX ) {
      float
	*buf_obj, *tempp, d_value1, d_value2;

      PLA_Obj_local_buffer( obj, (void **) &buf_obj );
      
      if ( datatype == MPI_DOUBLE ){
	d_value1 = ( double ) *( ( double *) value );
	d_value2 = 0.0;
      }
      else if ( datatype == MPI_FLOAT ){
      }
      else
	PLA_Abort( "Datatype not yet implemented", __LINE__, __FILE__ );

      for ( j=0; j<n; j++ ){
	tempp = buf_obj   + j*ldim*2;
	for ( i=0; i<m; i++ ){
	  *tempp++ = d_value1;
	  *tempp++ = d_value2;
	}
      }
    }
    else if ( datatype_obj == MPI_INT ) {
      int
	*buf_obj, *tempp, d_value;

      PLA_Obj_local_buffer( obj, (void **) &buf_obj );

      if ( datatype == MPI_DOUBLE )
	d_value = (int ) *( ( double *) value );
      else if ( datatype == MPI_FLOAT )
	d_value = ( int ) *( ( float *) value );
      else if ( datatype == MPI_INT )
	d_value = ( int ) *( ( int *) value );
      else
	PLA_Abort( "Datatype not yet implemented", __LINE__, __FILE__ );

      for ( j=0; j<n; j++ ){
	tempp = buf_obj   + j*ldim;
	for ( i=0; i<m; i++ )
	  *tempp++ = d_value;
      }
    }
    else {
      PLA_Abort( "Datatype not yet implemented", __LINE__, __FILE__ );
    }
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_set_exit( obj, datatype, value );

  return return_value;
}

/***************************************************************************/

int PLA_Obj_set_to_zero   (PLA_Obj      obj)

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to zero (0).
The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE, etc.)
where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS;

  double zero = 0.0;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_set_to_zero_enter( obj );

  PLA_Obj_set( obj, MPI_DOUBLE, &zero );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_set_to_zero_exit( obj );

  return return_value;
}

/***************************************************************************/

int PLA_Obj_set_to_one   (PLA_Obj      obj)

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to one (1).
The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE, etc.)
where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS;

  double one = 1.0;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_set_to_one_enter( obj );

  PLA_Obj_set( obj, MPI_DOUBLE, &one );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_set_to_one_exit( obj );

  return return_value;
}


/***************************************************************************/

int PLA_Obj_set_to_minus_one   (PLA_Obj      obj)

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to negative
one (-1). The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE,
etc.) where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS;

  double minus_one = -1.0;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_set_to_minus_one_enter( obj );

  PLA_Obj_set( obj, MPI_DOUBLE, &minus_one );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_set_to_minus_one_exit( obj );

  return return_value;
}


/* Chapter 3: Advanced Linear Algebra Object Manipulation

  Section 3.1: Creating Views     Section 3.5: Creating Objects
  into Objects                    "Conformal to" ...
  PLA_Obj_view                    PLA_Vector_create_conf_to
  PLA_Obj_view_all                PLA_Mvector_create_conf_to
  PLA_Obj_view_swap               PLA_Pvector_create_conf_to
                                  PLA_Pmvector_create_conf_to
                                  PLA_Matrix_create_conf_to
  Section 3.2: Splitting of       PLA_Mscalar_create_conf_to
  LAObjs                          PLA_Create_constants_conf_to
  PLA_Obj_split_4
  PLA_Obj_horz_split_2
  PLA_Obj_vert_split_2            Section 3.6: Annotating Object
                                  Orientation
                                  PLA_Obj_set_orientation
  Section 3.3: Shifting of        PLA_Obj_get_orientation
  LAObjs
  PLA_Obj_view_shift

                                  Section 3.7: Casting Object Types
  Section 3.4: Determining        PLA_Obj_objtype_cast
  Where to Split
  PLA_Obj_split_size

----------------------------------------------------------------------------*/

/***************************************************************************/

int PLA_Obj_view   (
                    PLA_Obj   obj,             int      global_length,
                    int         global_width,    int      align_row,
                    int         align_col,       PLA_Obj *new_obj)

/*----------------------------------------------------------------------------

Purpose : Create a view (reference) into an existing linear algebra
object.


IN     obj              object into which view is taken
IN     global_length    row dimension of view
IN     global_width     column dimension of view
IN     align_row        row index in old object of upper-left-hand element of view
IN     align_col        column index in old object of upper-left-hand element of view
IN/OUT new_obj          created view

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS,
    objtype;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_view_enter( 
         obj, global_length, global_width, align_row, align_col, new_obj );

  PLA_Obj_split_4( obj, 
		   ( align_row == PLA_ALIGN_FIRST ? 0 : align_row ), 
		   ( align_col == PLA_ALIGN_FIRST ? 0 : align_col ), 
		    PLA_DUMMY,   PLA_DUMMY,
		    PLA_DUMMY,   new_obj );

  if ( global_length == PLA_DIM_ALL )
    PLA_OBJ_GLOBAL_LENGTH( *new_obj, &global_length );
  if ( global_width == PLA_DIM_ALL )
    PLA_OBJ_GLOBAL_WIDTH( *new_obj, &global_width );

  PLA_Obj_split_4( *new_obj, global_length, global_width, 
		                   new_obj,     PLA_DUMMY,
		                   PLA_DUMMY,   PLA_DUMMY );
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_view_exit( 
         obj, global_length, global_width, align_row, align_col, new_obj );

  return return_value;
}

/***************************************************************************/

int PLA_Obj_view_all   (PLA_Obj     old_obj,      PLA_Obj     *new_obj)

/*----------------------------------------------------------------------------

Purpose : Create a view (reference) into all of existing linear algebra
object.


IN          old_obj          object into which view is taken
IN/OUT      new_obj          created view

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS,
    objtype,
    i;
  PLA_Obj
    view = NULL;
  char 
    *tempp_new, *tempp_old;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_view_all_enter( old_obj, new_obj );

  if ( old_obj != *new_obj ){
    if ( *new_obj != NULL )
      PLA_Obj_free( new_obj );

    view     = pla_get_view( );
  
    tempp_new = ( char * ) view;
    tempp_old = ( char * ) old_obj;

    memcpy( tempp_new, tempp_old, sizeof( struct PLA_Obj_view_struct ) );
/*    for ( i=0; i<sizeof( struct PLA_Obj_view_struct ); i++ )
      *tempp_new++ = *tempp_old++; */

    view->base_obj->number_of_views++;

    *new_obj = view;
  }

  PLA_Obj_objtype( old_obj, &objtype );
  if ( objtype == PLA_MATRIX ) {
    if ( view->proj_onto == PLA_PROJ_ONTO_ROW )
    view->proj_onto = PLA_PROJ_ONTO_COL;
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_view_all_exit( old_obj, new_obj );

  return return_value;
}


/***************************************************************************/

int PLA_Obj_view_swap    ( PLA_Obj     *obj1,        PLA_Obj   *obj2)

/*----------------------------------------------------------------------------

Purpose : Swaps two views.


IN/OUT      obj1          linear algebra object
IN/OUT      obj2          linear algebra object

----------------------------------------------------------------------------*/
{
  PLA_Obj 
    temp = NULL;

  temp = *obj1;
  *obj1 = *obj2;
  *obj2 = temp;

  return PLA_SUCCESS;
}

/***************************************************************************/

int PLA_Obj_split_4  (
     PLA_Obj       obj,
     int            length,            int           width,
     PLA_Obj       *upper_left_obj,   PLA_Obj       *upper_right_obj,
     PLA_Obj       *lower_left_obj,   PLA_Obj       *lower_right_obj)

/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into four quadrants

IN         obj               object to be split
IN         length            row dimension of block (as indicated below)
IN         width             column dimension of block (as indicated below)
OUT        upper_left_obj    new object for upper-left block
OUT        upper_right_obj   new object for upper-right block
OUT        lower_left_obj    new object for lower-left block
OUT        lower_right_obj   new object for lower-right block

Here the value of length and width determine the size of one block while
their sign determines which block the length and widht refer to.
If :
length >= 0        width >= 0            upper_left_obj
length <  0        width >= 0            lower_left_obj
length >= 0        width <  0            upper_right_obj
length <  0        width <  0            lower_right_obj

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS,
    objtype;


  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_split_4_enter( 
	  obj, length, width,  upper_left_obj, upper_right_obj,
                               lower_left_obj, lower_right_obj );

  PLA_OBJ_OBJTYPE( obj, &objtype );

  switch( objtype ){
  case PLA_MVECTOR:
    PLA_Obj_split_4_mv( 
	  obj, length, width,  upper_left_obj, upper_right_obj,
                               lower_left_obj, lower_right_obj );
    break;
  case PLA_MATRIX:
    PLA_Obj_split_4_matrix(
	  obj, length, width,  upper_left_obj, upper_right_obj,
                               lower_left_obj, lower_right_obj ); 
    break;
  case PLA_MSCALAR:
    PLA_Obj_split_4_mscalar(
	  obj, length, width,  upper_left_obj, upper_right_obj,
                               lower_left_obj, lower_right_obj );
    break;
  case PLA_PMVECTOR:
    PLA_Obj_split_4_pmvector(
	  obj, length, width,  upper_left_obj, upper_right_obj,
                               lower_left_obj, lower_right_obj ); 
    break; 
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_split_4_exit( 
	  obj, length, width,  upper_left_obj, upper_right_obj,
                               lower_left_obj, lower_right_obj );

  return return_value;
}

/***************************************************************************/

int PLA_Obj_split_4_mv  (
     PLA_Obj       obj,
     int            length,            int           width,
     PLA_Obj       *upper_left_obj,   PLA_Obj       *upper_right_obj,
     PLA_Obj       *lower_left_obj,   PLA_Obj       *lower_right_obj)

/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into four quadrants

IN         obj               object to be split
IN         length            row dimension of block (as indicated below)
IN         width             column dimension of block (as indicated below)
OUT        upper_left_obj    new object for upper-left block
OUT        upper_right_obj   new object for upper-right block
OUT        lower_left_obj    new object for lower-left block
OUT        lower_right_obj   new object for lower-right block

Here the value of length and width determine the size of one block while
their sign determines which block the length and widht refer to.
If :
length >= 0        width >= 0            upper_left_obj
length <  0        width >= 0            lower_left_obj
length >= 0        width <  0            upper_right_obj
length <  0        width <  0            lower_right_obj

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS,
    global_length, upper_global_length, lower_global_length,
    global_width,  left_global_width,   right_global_width,
    global_align_row, upper_global_align_row, lower_global_align_row,
    upper_local_length, lower_local_length,
    left_local_width,   left_local_length,
    upper_size_top, lower_size_top, upper_size_bottom, lower_size_bottom,
    upper_owner_top, lower_owner_top, upper_owner_bottom, lower_owner_bottom,
    typesize, local_ldim;

  MPI_Datatype 
    datatype;

  PLA_Obj
    view = NULL; 
  
  PLA_Template 
    templ = NULL;

  if ( upper_left_obj != PLA_DUMMY && *upper_left_obj != obj )
    PLA_Obj_view_all( obj, upper_left_obj );

  if ( upper_right_obj != PLA_DUMMY && *upper_right_obj != obj )
    PLA_Obj_view_all( obj, upper_right_obj );

  if ( lower_left_obj != PLA_DUMMY && *lower_left_obj != obj )
    PLA_Obj_view_all( obj, lower_left_obj );

  if ( lower_right_obj != PLA_DUMMY && *lower_right_obj != obj )
    PLA_Obj_view_all( obj, lower_right_obj );

  PLA_OBJ_GLOBAL_LENGTH   ( obj, &global_length );
  PLA_OBJ_GLOBAL_WIDTH    ( obj, &global_width );
  PLA_OBJ_GLOBAL_ALIGN_ROW( obj, &global_align_row );
  PLA_OBJ_TEMPLATE        ( obj, &templ );
  PLA_OBJ_LOCAL_LDIM      ( obj, &local_ldim );
  PLA_OBJ_DATATYPE        ( obj, &datatype );

  MPI_Type_size( datatype, &typesize );

  if ( length >= 0 ){
    upper_global_length     = length;
    lower_global_length     = global_length - length;
    upper_global_align_row  = global_align_row;
    lower_global_align_row  = global_align_row + length;
  }
  else {
    upper_global_length     = global_length + length;
    lower_global_length     = -length;
    upper_global_align_row  = global_align_row;
    lower_global_align_row  = global_align_row + upper_global_length;
  }

  if ( width >= 0 ){
    left_global_width       = width;
    right_global_width      = global_width - width;
  }
  else {
    left_global_width       = global_width + width;
    right_global_width      = -width;
  }

  PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ALL,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       upper_global_align_row,        /* Alignment to template wrt direction */
       upper_global_length,           /* Length in direction */ 
       &upper_local_length,           /* Length assigned to this node */
       &upper_size_top,               /* Size of first block */
       &upper_owner_top,             /* Owner of first block */
       &upper_size_bottom,             /* Size of last block */
       &upper_owner_bottom);           /* Owner of last block */

  PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ALL,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       lower_global_align_row,       /* Alignment to template wrt direction */
       lower_global_length,          /* Length in direction */ 
       &lower_local_length,          /* Length assigned to this node */
       &lower_size_top,              /* Size of first block */
       &lower_owner_top,	      /* Owner of first block */
       &lower_size_bottom,           /* Size of last block */
       &lower_owner_bottom );        /* Owner of last block */

  if ( upper_left_obj != PLA_DUMMY ){
    view = *upper_left_obj;
    view->global_length     = upper_global_length;
    view->global_width      = left_global_width;
    view->global_align_row  = upper_global_align_row;

    view->local_length      = upper_local_length;
    view->local_width       = left_global_width;

    view->size_top          = upper_size_top;
    view->size_bottom       = upper_size_bottom;
    view->size_left         = left_global_width;
    view->size_right        = left_global_width;
    
 /* view->local_buffer       = unchanged */
  }

  if ( upper_right_obj != PLA_DUMMY ){
    view = *upper_right_obj;
    view->global_length     = upper_global_length;
    view->global_width      = right_global_width;
    view->global_align_row  = upper_global_align_row;

    view->global_align_col += width;   /* Consistent with R12, but fix!!! */

    view->local_length      = upper_local_length;
    view->local_width       = right_global_width;

    view->size_top          = upper_size_top;
    view->owner_top         = upper_owner_top;
    view->size_bottom       = upper_size_bottom;
    view->owner_bottom      = upper_owner_bottom;
    view->size_left         = right_global_width;
    view->size_right        = right_global_width;
    
    view->local_buffer     = (char *) view->local_buffer +
                              left_global_width * local_ldim * typesize;
  }

  if ( lower_left_obj != PLA_DUMMY ){
    view = *lower_left_obj;
    view->global_length     = lower_global_length;
    view->global_width      = left_global_width;
    view->global_align_row  = lower_global_align_row;

    view->local_length      = lower_local_length;
    view->local_width       = left_global_width;

    view->size_top          = lower_size_top;
    view->owner_top         = lower_owner_top;
    view->size_bottom       = lower_size_bottom;
    view->owner_bottom      = lower_owner_bottom;
    view->size_left         = left_global_width;
    view->size_right        = left_global_width;
    
    view->local_buffer     =  (char *) view->local_buffer + 
                              upper_local_length * typesize;
  }

  if ( lower_right_obj != PLA_DUMMY ){
    view = *lower_right_obj;
    view->global_length     = lower_global_length;
    view->global_width      = right_global_width;
    view->global_align_row  = lower_global_align_row;

    view->global_align_col += width;   /* Consistent with R12, but fix!!! */

    view->local_length      = lower_local_length;
    view->local_width       = right_global_width;

    view->size_top          = lower_size_top;
    view->owner_top         = lower_owner_top;
    view->size_bottom       = lower_size_bottom;
    view->owner_bottom      = lower_owner_bottom;
    view->size_left         = right_global_width;
    view->size_right        = right_global_width;
    
    view->local_buffer     =  (char *) view->local_buffer + 
                              left_global_width  * local_ldim * typesize +
                              upper_local_length * typesize;
  }


  return return_value;
}

/***************************************************************************/

int PLA_Obj_split_4_matrix  (
     PLA_Obj       obj,
     int            length,            int           width,
     PLA_Obj       *upper_left_obj,   PLA_Obj       *upper_right_obj,
     PLA_Obj       *lower_left_obj,   PLA_Obj       *lower_right_obj)

/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into four quadrants

IN         obj               object to be split
IN         length            row dimension of block (as indicated below)
IN         width             column dimension of block (as indicated below)
OUT        upper_left_obj    new object for upper-left block
OUT        upper_right_obj   new object for upper-right block
OUT        lower_left_obj    new object for lower-left block
OUT        lower_right_obj   new object for lower-right block

Here the value of length and width determine the size of one block while
their sign determines which block the length and widht refer to.
If :
length >= 0        width >= 0            upper_left_obj
length <  0        width >= 0            lower_left_obj
length >= 0        width <  0            upper_right_obj
length <  0        width <  0            lower_right_obj

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS,
    global_length, upper_global_length, lower_global_length,
    global_width,  left_global_width,   right_global_width,
    upper_local_length, lower_local_length,
    left_local_width,   right_local_width,
    global_align_row, upper_global_align_row, lower_global_align_row,
    global_align_col, left_global_align_col, right_global_align_col,
    upper_size_top, lower_size_top, upper_size_bottom, lower_size_bottom,
    upper_owner_top, lower_owner_top, upper_owner_bottom, lower_owner_bottom,
    left_size_left, right_size_left, left_size_right, right_size_right,
    left_owner_left, right_owner_left, left_owner_right, right_owner_right,
    typesize, local_ldim;

  MPI_Datatype 
    datatype;

  PLA_Obj
    view = NULL; 
  
  PLA_Template 
    templ = NULL;

  if ( upper_left_obj != PLA_DUMMY && *upper_left_obj != obj )
    PLA_Obj_view_all( obj, upper_left_obj );

  if ( upper_right_obj != PLA_DUMMY && *upper_right_obj != obj )
    PLA_Obj_view_all( obj, upper_right_obj );

  if ( lower_left_obj != PLA_DUMMY && *lower_left_obj != obj )
    PLA_Obj_view_all( obj, lower_left_obj );

  if ( lower_right_obj != PLA_DUMMY && *lower_right_obj != obj )
    PLA_Obj_view_all( obj, lower_right_obj );

  PLA_OBJ_GLOBAL_LENGTH   ( obj, &global_length );
  PLA_OBJ_GLOBAL_WIDTH    ( obj, &global_width );
  PLA_OBJ_GLOBAL_ALIGN_ROW( obj, &global_align_row );
  PLA_OBJ_GLOBAL_ALIGN_COL( obj, &global_align_col );
  PLA_OBJ_TEMPLATE        ( obj, &templ );
  PLA_OBJ_LOCAL_LDIM      ( obj, &local_ldim );
  PLA_OBJ_DATATYPE        ( obj, &datatype );

  MPI_Type_size( datatype, &typesize );

  if ( length >= 0 ){
    upper_global_length     = length;
    lower_global_length     = global_length - length;
    upper_global_align_row  = global_align_row;
    lower_global_align_row  = global_align_row + length;
  }
  else {
    upper_global_length     = global_length + length;
    lower_global_length     = -length;
    upper_global_align_row  = global_align_row;
    lower_global_align_row  = global_align_row + upper_global_length;
  }

  if ( width >= 0 ){
    left_global_width       = width;
    right_global_width      = global_width - width;
    left_global_align_col   = global_align_col;
    right_global_align_col  = global_align_col + width;
  }
  else {
    left_global_width       = global_width + width;
    right_global_width      = -width;
    left_global_align_col   = global_align_col;
    right_global_align_col  = global_align_col + left_global_width;
  }

  PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_COL,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       upper_global_align_row,        /* Alignment to template wrt direction */
       upper_global_length,           /* Length in direction */ 
       &upper_local_length,           /* Length assigned to this node */
       &upper_size_top,               /* Size of first block */
       &upper_owner_top,             /* Owner of first block */
       &upper_size_bottom,             /* Size of last block */
       &upper_owner_bottom);           /* Owner of last block */

  PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_COL,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       lower_global_align_row,       /* Alignment to template wrt direction */
       lower_global_length,          /* Length in direction */ 
       &lower_local_length,          /* Length assigned to this node */
       &lower_size_top,              /* Size of first block */
       &lower_owner_top,	      /* Owner of first block */
       &lower_size_bottom,           /* Size of last block */
       &lower_owner_bottom );        /* Owner of last block */

  PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ROW,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       left_global_align_col,         /* Alignment to template wrt direction */
       left_global_width,             /* Length in direction */ 
       &left_local_width,             /* Length assigned to this node */
       &left_size_left,               /* Size of first block */
       &left_owner_left,              /* Owner of first block */
       &left_size_right,              /* Size of last block */
       &left_owner_right);            /* Owner of last block */

  PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ROW,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       right_global_align_col,        /* Alignment to template wrt direction */
       right_global_width,            /* Length in direction */ 
       &right_local_width,            /* Length assigned to this node */
       &right_size_left,              /* Size of first block */
       &right_owner_left,	      /* Owner of first block */
       &right_size_right,             /* Size of last block */
       &right_owner_right );          /* Owner of last block */

  if ( upper_left_obj != PLA_DUMMY ){
    view = *upper_left_obj;
    view->global_length     = upper_global_length;
    view->global_width      = left_global_width;
    view->global_align_row  = upper_global_align_row;
    view->global_align_col  = left_global_align_col;

    view->local_length      = upper_local_length;
    view->local_width       = left_local_width;

    view->size_top          = upper_size_top;
    view->owner_top         = upper_owner_top;
    view->size_bottom       = upper_size_bottom;
    view->owner_bottom      = upper_owner_bottom;
    view->size_left         = left_size_left;
    view->owner_left        = left_owner_left;
    view->size_right        = left_size_right;
    view->owner_right       = left_owner_right;
 /* view->local_buffer       = unchanged */
    view->proj_onto         = PLA_PROJ_ONTO_COL;
  }

  if ( upper_right_obj != PLA_DUMMY ){
    view = *upper_right_obj;
    view->global_length     = upper_global_length;
    view->global_width      = right_global_width;
    view->global_align_row  = upper_global_align_row;
    view->global_align_col  = right_global_align_col;

    view->local_length      = upper_local_length;
    view->local_width       = right_local_width;

    view->size_top          = upper_size_top;
    view->owner_top         = upper_owner_top;
    view->size_bottom       = upper_size_bottom;
    view->owner_bottom      = upper_owner_bottom;
    view->size_left         = right_size_left;
    view->owner_left        = right_owner_left;
    view->size_right        = right_size_right;
    view->owner_right       = right_owner_right;
    
    view->local_buffer     = (char *) view->local_buffer +
                              left_local_width * local_ldim * typesize;
    view->proj_onto         = PLA_PROJ_ONTO_COL;
  }

  if ( lower_left_obj != PLA_DUMMY ){
    view = *lower_left_obj;
    view->global_length     = lower_global_length;
    view->global_width      = left_global_width;
    view->global_align_row  = lower_global_align_row;
    view->global_align_col  = left_global_align_col;

    view->local_length      = lower_local_length;
    view->local_width       = left_local_width;

    view->size_top          = lower_size_top;
    view->owner_top         = lower_owner_top;
    view->size_bottom       = lower_size_bottom;
    view->owner_bottom      = lower_owner_bottom;
    view->size_left         = left_size_left;
    view->owner_left        = left_owner_left;
    view->size_right        = left_size_right;
    view->owner_right       = left_owner_right;
    
    view->local_buffer     =  (char *) view->local_buffer + 
                              upper_local_length * typesize;
    view->proj_onto         = PLA_PROJ_ONTO_COL;
  }

  if ( lower_right_obj != PLA_DUMMY ){
    view = *lower_right_obj;
    view->global_length     = lower_global_length;
    view->global_width      = right_global_width;
    view->global_align_row  = lower_global_align_row;
    view->global_align_col  = right_global_align_col;

    view->local_length      = lower_local_length;
    view->local_width       = right_local_width;

    view->size_top          = lower_size_top;
    view->owner_top         = lower_owner_top;
    view->size_bottom       = lower_size_bottom;
    view->owner_bottom      = lower_owner_bottom;
    view->size_left         = right_size_left;
    view->owner_left        = right_owner_left;
    view->size_right        = right_size_right;
    view->owner_right       = right_owner_right;
    
    view->local_buffer      =  (char *) view->local_buffer + 
                               left_local_width  * local_ldim * typesize +
                               upper_local_length * typesize;
    view->proj_onto         = PLA_PROJ_ONTO_COL;
  }

  return return_value;
}

/***************************************************************************/

int PLA_Obj_split_4_mscalar  (
     PLA_Obj       obj,
     int            length,            int           width,
     PLA_Obj       *upper_left_obj,   PLA_Obj       *upper_right_obj,
     PLA_Obj       *lower_left_obj,   PLA_Obj       *lower_right_obj)

/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into four quadrants

IN         obj               object to be split
IN         length            row dimension of block (as indicated below)
IN         width             column dimension of block (as indicated below)
OUT        upper_left_obj    new object for upper-left block
OUT        upper_right_obj   new object for upper-right block
OUT        lower_left_obj    new object for lower-left block
OUT        lower_right_obj   new object for lower-right block

Here the value of length and width determine the size of one block while
their sign determines which block the length and widht refer to.
If :
length >= 0        width >= 0            upper_left_obj
length <  0        width >= 0            lower_left_obj
length >= 0        width <  0            upper_right_obj
length <  0        width <  0            lower_right_obj

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS,
    global_length, upper_global_length, lower_global_length,
    global_width,  left_global_width,   right_global_width,
    typesize, local_ldim,
    myrow, mycol, owner_row, owner_col;

  MPI_Datatype 
    datatype;

  PLA_Obj
    view = NULL; 

  PLA_Template
    templ = NULL;
  
  if ( upper_left_obj != PLA_DUMMY && *upper_left_obj != obj )
    PLA_Obj_view_all( obj, upper_left_obj );

  if ( upper_right_obj != PLA_DUMMY && *upper_right_obj != obj )
    PLA_Obj_view_all( obj, upper_right_obj );

  if ( lower_left_obj != PLA_DUMMY && *lower_left_obj != obj )
    PLA_Obj_view_all( obj, lower_left_obj );

  if ( lower_right_obj != PLA_DUMMY && *lower_right_obj != obj )
    PLA_Obj_view_all( obj, lower_right_obj );

  PLA_OBJ_GLOBAL_LENGTH   ( obj, &global_length );
  PLA_OBJ_GLOBAL_WIDTH    ( obj, &global_width );
  PLA_OBJ_LOCAL_LDIM      ( obj, &local_ldim );
  PLA_OBJ_DATATYPE        ( obj, &datatype );

  PLA_OBJ_TEMPLATE        ( obj, &templ );
  PLA_Temp_comm_col_rank   ( templ, &myrow );
  PLA_Temp_comm_row_rank   ( templ, &mycol );

  PLA_OBJ_OWNER_ROW       ( obj, &owner_row );
  PLA_OBJ_OWNER_COL       ( obj, &owner_col );

  MPI_Type_size( datatype, &typesize );

  if ( length >= 0 ){
    upper_global_length     = length;
    lower_global_length     = global_length - length;
  }
  else {
    upper_global_length     = global_length + length;
    lower_global_length     = -length;
  }

  if ( width >= 0 ){
    left_global_width       = width;
    right_global_width      = global_width - width;
  }
  else {
    left_global_width       = global_width + width;
    right_global_width      = -width;
  }

  if ( upper_left_obj != PLA_DUMMY ){
    view = *upper_left_obj;
    view->global_length     = upper_global_length;
    view->global_width      = left_global_width;

    if ( ( owner_row == myrow || owner_row == PLA_ALL_ROWS ) &&
         ( owner_col == mycol || owner_col == PLA_ALL_COLS ) ) {
      view->local_length      = upper_global_length;
      view->local_width       = left_global_width;
    }

    view->size_top          = upper_global_length;
 /* view->owner_top         = unchanged */
    view->size_bottom       = upper_global_length;
 /* view->owner_bottom      = unchanged */
    view->size_left         = left_global_width;
 /* view->owner_left        = unchanged */
    view->size_right        = left_global_width;
 /* view->owner_right       = unchanged */
 /* view->local_buffer       = unchanged */
  }

  if ( upper_right_obj != PLA_DUMMY ){
    view = *upper_right_obj;
    view->global_length     = upper_global_length;
    view->global_width      = right_global_width;
    view->global_align_col += left_global_width;   /* consistent w/ R12 */

    if ( ( owner_row == myrow || owner_row == PLA_ALL_ROWS ) &&
         ( owner_col == mycol || owner_col == PLA_ALL_COLS ) ) {
      view->local_length      = upper_global_length;
      view->local_width       = right_global_width;
      view->local_buffer     = (char *) view->local_buffer +
                               left_global_width * local_ldim * typesize;
    }

    view->size_top          = upper_global_length;
 /* view->owner_top         = unchanged */
    view->size_bottom       = upper_global_length;
 /* view->owner_bottom      = unchanged */
    view->size_left         = right_global_width;
 /* view->owner_left        = unchanged */
    view->size_right        = right_global_width;
 /* view->owner_right       = unchanged */
    
  }

  if ( lower_left_obj != PLA_DUMMY ){
    view = *lower_left_obj;
    view->global_length     = lower_global_length;
    view->global_width      = left_global_width;
    view->global_align_row += upper_global_length;  /* consistent w/ R12 */

    if ( ( owner_row == myrow || owner_row == PLA_ALL_ROWS ) &&
         ( owner_col == mycol || owner_col == PLA_ALL_COLS ) ) {
      view->local_length      = lower_global_length;
      view->local_width       = left_global_width;

      view->local_buffer     =  (char *) view->local_buffer + 
                                upper_global_length * typesize;
    }

    view->size_top          = lower_global_length;
 /* view->owner_top         = unchanged */
    view->size_bottom       = lower_global_length;
 /* view->owner_bottom      = unchanged */
    view->size_left         = left_global_width;
 /* view->owner_left        = unchanged */
    view->size_right        = left_global_width;
 /* view->owner_right       = unchanged */
  }

  if ( lower_right_obj != PLA_DUMMY ){
    view = *lower_right_obj;
    view->global_length     = lower_global_length;
    view->global_width      = right_global_width;
    view->global_align_row += upper_global_length; /* consistent w/ R12 */
    view->global_align_col += left_global_width;   /* consistent w/ R12 */

    if ( ( owner_row == myrow || owner_row == PLA_ALL_ROWS ) &&
         ( owner_col == mycol || owner_col == PLA_ALL_COLS ) ) {
      view->local_length      = lower_global_length;
      view->local_width       = right_global_width;
      view->local_buffer      =  (char *) view->local_buffer + 
                                 left_global_width  * local_ldim * typesize +
                                 upper_global_length * typesize;
    }

    view->size_top          = lower_global_length;
 /* view->owner_top         = unchanged */
    view->size_bottom       = lower_global_length;
 /* view->owner_bottom      = unchanged */
    view->size_left         = right_global_width;
 /* view->owner_left        = unchanged */
    view->size_right        = right_global_width;
 /* view->owner_right       = unchanged */
   
  }

  return return_value;
}



int PLA_Obj_split_4_pmvector( 
     PLA_Obj       obj,
     int            length,            int           width,
     PLA_Obj       *upper_left_obj,   PLA_Obj       *upper_right_obj,
     PLA_Obj       *lower_left_obj,   PLA_Obj       *lower_right_obj)

/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into four quadrants

IN         obj               object to be split
IN         length            row dimension of block (as indicated below)
IN         width             column dimension of block (as indicated below)
OUT        upper_left_obj    new object for upper-left block
OUT        upper_right_obj   new object for upper-right block
OUT        lower_left_obj    new object for lower-left block
OUT        lower_right_obj   new object for lower-right block

Here the value of length and width determine the size of one block while
their sign determines which block the length and widht refer to.
If :
length >= 0        width >= 0            upper_left_obj
length <  0        width >= 0            lower_left_obj
length >= 0        width <  0            upper_right_obj
length <  0        width <  0            lower_right_obj

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS,
    global_length, upper_global_length, lower_global_length,
    global_width,  left_global_width,   right_global_width,
    upper_local_length, lower_local_length,
    left_local_width,   right_local_width,
    global_align_row, upper_global_align_row, lower_global_align_row,
    global_align_col, left_global_align_col, right_global_align_col,
    upper_size_top, lower_size_top, upper_size_bottom, lower_size_bottom,
    upper_owner_top, lower_owner_top, upper_owner_bottom, lower_owner_bottom,
    left_size_left, right_size_left, left_size_right, right_size_right,
    left_owner_left, right_owner_left, left_owner_right, right_owner_right,
    typesize, local_ldim,
    myrow, mycol, owner_row, owner_col, project_onto;

  MPI_Datatype 
    datatype;

  PLA_Obj
    view = NULL; 
  
  PLA_Template 
    templ = NULL;

  if ( upper_left_obj != PLA_DUMMY && *upper_left_obj != obj )
    PLA_Obj_view_all( obj, upper_left_obj );

  if ( upper_right_obj != PLA_DUMMY && *upper_right_obj != obj )
    PLA_Obj_view_all( obj, upper_right_obj );

  if ( lower_left_obj != PLA_DUMMY && *lower_left_obj != obj )
    PLA_Obj_view_all( obj, lower_left_obj );

  if ( lower_right_obj != PLA_DUMMY && *lower_right_obj != obj )
    PLA_Obj_view_all( obj, lower_right_obj );

  PLA_OBJ_GLOBAL_LENGTH   ( obj, &global_length );
  PLA_OBJ_GLOBAL_WIDTH    ( obj, &global_width );
  PLA_OBJ_GLOBAL_ALIGN_ROW( obj, &global_align_row );
  PLA_OBJ_GLOBAL_ALIGN_COL( obj, &global_align_col );
  PLA_OBJ_OWNER_ROW       ( obj, &owner_row );
  PLA_OBJ_OWNER_COL       ( obj, &owner_col );
  PLA_OBJ_PROJECT_ONTO    ( obj, &project_onto );
  PLA_OBJ_TEMPLATE        ( obj, &templ );
  PLA_OBJ_LOCAL_LDIM      ( obj, &local_ldim );
  PLA_OBJ_DATATYPE        ( obj, &datatype );

  PLA_Temp_comm_col_rank( templ, &myrow );
  PLA_Temp_comm_row_rank( templ, &mycol );

  MPI_Type_size( datatype, &typesize );

  if ( length >= 0 ){
    upper_global_length     = length;
    lower_global_length     = global_length - length;
    upper_global_align_row  = global_align_row;
    lower_global_align_row  = global_align_row + length;
  }
  else {
    upper_global_length     = global_length + length;
    lower_global_length     = -length;
    upper_global_align_row  = global_align_row;
    lower_global_align_row  = global_align_row + upper_global_length;
  }

  if ( width >= 0 ){
    left_global_width       = width;
    right_global_width      = global_width - width;
    left_global_align_col   = global_align_col;
    right_global_align_col  = global_align_col + width;
  }
  else {
    left_global_width       = global_width + width;
    right_global_width      = -width;
    left_global_align_col   = global_align_col;
    right_global_align_col  = global_align_col + left_global_width;
  }

  if ( project_onto == PLA_PROJ_ONTO_COL ) {
    PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_COL,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       upper_global_align_row,        /* Alignment to template wrt direction */
       upper_global_length,           /* Length in direction */ 
       &upper_local_length,           /* Length assigned to this node */
       &upper_size_top,               /* Size of first block */
       &upper_owner_top,             /* Owner of first block */
       &upper_size_bottom,             /* Size of last block */
       &upper_owner_bottom);           /* Owner of last block */

    PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_COL,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       lower_global_align_row,       /* Alignment to template wrt direction */
       lower_global_length,          /* Length in direction */ 
       &lower_local_length,          /* Length assigned to this node */
       &lower_size_top,              /* Size of first block */
       &lower_owner_top,	      /* Owner of first block */
       &lower_size_bottom,           /* Size of last block */
       &lower_owner_bottom );        /* Owner of last block */
  }
  else {
    PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ROW,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       left_global_align_col,         /* Alignment to template wrt direction */
       left_global_width,             /* Length in direction */ 
       &left_local_width,             /* Length assigned to this node */
       &left_size_left,               /* Size of first block */
       &left_owner_left,              /* Owner of first block */
       &left_size_right,              /* Size of last block */
       &left_owner_right);            /* Owner of last block */

    PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ROW,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       right_global_align_col,        /* Alignment to template wrt direction */
       right_global_width,            /* Length in direction */ 
       &right_local_width,            /* Length assigned to this node */
       &right_size_left,              /* Size of first block */
       &right_owner_left,	      /* Owner of first block */
       &right_size_right,             /* Size of last block */
       &right_owner_right );          /* Owner of last block */
  }

  if ( upper_left_obj != PLA_DUMMY ){
    view = *upper_left_obj;
    view->global_length     = upper_global_length;
    view->global_width      = left_global_width;
    view->global_align_row  = upper_global_align_row;
    view->global_align_col  = left_global_align_col;

    if ( project_onto == PLA_PROJ_ONTO_COL ){
      view->local_length      = upper_local_length;

      if ( mycol == owner_col || PLA_ALL_COLS == owner_col ) {
	view->local_width       = left_global_width;
      }

      view->size_top          = upper_size_top;
      view->owner_top         = upper_owner_top;
      view->size_bottom       = upper_size_bottom;
      view->owner_bottom      = upper_owner_bottom;
      view->size_left         = left_global_width;
/*      view->owner_left        = left_owner_left; */
      view->size_right        = left_global_width;
/*      view->owner_right       = left_owner_right; */

      /* view->local_buffer       = unchanged */
    }
    else {
      if ( myrow == owner_row || PLA_ALL_ROWS == owner_row ) {
	view->local_length      = upper_global_length;
      }

      view->local_width       = left_local_width;

      view->size_top          = upper_global_length;
/*      view->owner_top         = upper_owner_top; */
      view->size_bottom       = upper_global_length;
/*      view->owner_bottom      = upper_owner_bottom; */
      view->size_left         = left_size_left;
      view->owner_left        = left_owner_left; 
      view->size_right        = left_size_right;
      view->owner_right       = left_owner_right; 

      /* view->local_buffer       = unchanged */
    }
  }

  if ( upper_right_obj != PLA_DUMMY ){
    view = *upper_right_obj;
    view->global_length     = upper_global_length;
    view->global_width      = right_global_width;
    view->global_align_row  = upper_global_align_row;
    view->global_align_col  = right_global_align_col;

    if ( project_onto == PLA_PROJ_ONTO_COL ){
      view->local_length      = upper_local_length;

      if ( mycol == owner_col || PLA_ALL_COLS == owner_col ) {
	view->local_width       = right_global_width;

	view->local_buffer     = (char *) view->local_buffer +
                  	left_global_width * local_ldim * typesize;
      }

      view->size_top          = upper_size_top;
      view->owner_top         = upper_owner_top;
      view->size_bottom       = upper_size_bottom;
      view->owner_bottom      = upper_owner_bottom;
      view->size_left         = right_global_width;
/*      view->owner_left        = right_owner_left; */
      view->size_right        = right_global_width;
/*      view->owner_right       = right_owner_right; */
    }
    else {
      if ( myrow == owner_row || PLA_ALL_ROWS == owner_row ) {
	view->local_length    = upper_global_length;

	view->local_buffer     = (char *) view->local_buffer +
                  	left_local_width * local_ldim * typesize;
      }

      view->local_width       = right_local_width;

      view->size_top          = upper_global_length;
/*      view->owner_top         = upper_owner_top; */
      view->size_bottom       = upper_global_length;
/*      view->owner_bottom      = upper_owner_bottom; */
      view->size_left         = right_size_left;
      view->owner_left        = right_owner_left;
      view->size_right        = right_size_right;
      view->owner_right       = right_owner_right; 
    }
  }
 
/*  continue editing here! */
  if ( lower_left_obj != PLA_DUMMY ){
    view = *lower_left_obj;
    view->global_length     = lower_global_length;
    view->global_width      = left_global_width;
    view->global_align_row  = lower_global_align_row;
    view->global_align_col  = left_global_align_col;

    if ( project_onto == PLA_PROJ_ONTO_COL ){
      view->local_length      = lower_local_length;

      if ( mycol == owner_col || PLA_ALL_COLS == owner_col ) {
	view->local_width       = left_global_width;

	view->local_buffer     =  (char *) view->local_buffer + 
                                upper_local_length * typesize;
      }

      view->size_top          = lower_size_top;
      view->owner_top         = lower_owner_top;
      view->size_bottom       = lower_size_bottom;
      view->owner_bottom      = lower_owner_bottom;
      view->size_left         = left_global_width;
/*      view->owner_left        = left_owner_left; */
      view->size_right        = left_global_width;
/*      view->owner_right       = left_owner_right; */

      /* view->local_buffer       = unchanged */
    }
    else {
      if ( myrow == owner_row || PLA_ALL_ROWS == owner_row ) {
	view->local_length      = lower_global_length;

	view->local_buffer     =  (char *) view->local_buffer + 
                                upper_global_length * typesize;
      }

      view->local_width       = left_local_width;

      view->size_top          = lower_global_length;
/*      view->owner_top         = lower_owner_top; */
      view->size_bottom       = lower_global_length;
/*      view->owner_bottom      = lower_owner_bottom; */
      view->size_left         = left_size_left;
      view->owner_left        = left_owner_left; 
      view->size_right        = left_size_right;
      view->owner_right       = left_owner_right; 

    }
  }

  if ( lower_right_obj != PLA_DUMMY ){
    view = *lower_right_obj;
    view->global_length     = lower_global_length;
    view->global_width      = right_global_width;
    view->global_align_row  = lower_global_align_row;
    view->global_align_col  = right_global_align_col;

    if ( project_onto == PLA_PROJ_ONTO_COL ){
      view->local_length      = lower_local_length;

      if ( mycol == owner_col || PLA_ALL_COLS == owner_col ) {
	view->local_width       = right_global_width;

	view->local_buffer      =  (char *) view->local_buffer + 
                               left_global_width  * local_ldim * typesize +
                               upper_local_length * typesize;
      }

      view->size_top          = lower_size_top;
      view->owner_top         = lower_owner_top;
      view->size_bottom       = lower_size_bottom;
      view->owner_bottom      = lower_owner_bottom;
      view->size_left         = right_global_width;
/*      view->owner_left        = right_owner_left; */
      view->size_right        = right_global_width;
/*      view->owner_right       = right_owner_right; */
    }
    else {
      if ( myrow == owner_row || PLA_ALL_ROWS == owner_row ) {
	view->local_length    = lower_global_length;

	view->local_buffer      =  (char *) view->local_buffer + 
                               left_local_width  * local_ldim * typesize +
                               upper_global_length * typesize;
      }

      view->local_width       = right_local_width;

      view->size_top          = lower_global_length;
/*      view->owner_top         = lower_owner_top; */
      view->size_bottom       = lower_global_length;
/*      view->owner_bottom      = lower_owner_bottom; */
      view->size_left         = right_size_left;
      view->owner_left        = right_owner_left;
      view->size_right        = right_size_right;
      view->owner_right       = right_owner_right; 
    }
  }

  return return_value;
}

/***************************************************************************/

int PLA_Obj_horz_split_2   (  PLA_Obj     obj,      int          length,
                                                    PLA_Obj      *upper_obj,
                                                    PLA_Obj      *lower_obj)

/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into top and bottom


IN         obj               object to be split
IN         length            row dimension of block (as indicated below)
OUT        upper_obj         new object for upper block
OUT        lower_obj         new object for lower block

If :
length >= 0     length specifies   upper_obj
length <  0     length specifies   lower_obj

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    width;
  
  PLA_OBJ_GLOBAL_WIDTH( obj, &width );

  value = PLA_Obj_split_4( obj, length, width, upper_obj, PLA_DUMMY,
			                       lower_obj, PLA_DUMMY );

  return value;
}

/***************************************************************************/

int PLA_Obj_vert_split_2   (  PLA_Obj     obj,           int          width,
                              PLA_Obj      *left_obj,    PLA_Obj      *right_obj)

/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into left and right


IN         obj               object to be split
IN         width             column dimension of block (as indicated below)
OUT        left_obj          new object for left block
OUT        right_obj         new object for right block

If :
width >= 0     width specifies   left_obj
width <  0     width specifies   right_obj

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    length;

  PLA_OBJ_GLOBAL_LENGTH( obj, &length );

  value = PLA_Obj_split_4( obj, length, width, left_obj,  right_obj,
		  	                        PLA_DUMMY, PLA_DUMMY );

  return value;
}

/***************************************************************************/

int PLA_Obj_view_shift  (PLA_Obj      obj,           int      length_top,
                          int          width_left,    int      width_right,
                          int          length_bottom)

/*----------------------------------------------------------------------------

Purpose : Shift the boundaries of the view of the linear algebra object



IN/OUT       obj           object to be grown [shifted]
IN           length_top    length to shift top boundary
IN           width_left    width to shift left boundary
IN           width_right   width to shift right boundary
IN           length_bottom length to shift bottom boundary

Note : Positive is down/right negative is up/left (as if you were in quadrant IV
in the cartesian plane.

----------------------------------------------------------------------------*/
{
  int 
    return_value = PLA_SUCCESS,
    objtype;


  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_view_shift_enter( 
	  obj, length_top, width_left, width_right, length_bottom );

  PLA_OBJ_OBJTYPE( obj, &objtype );

  switch( objtype ){
  case PLA_MVECTOR:
    PLA_Obj_view_shift_mv( 
	  obj, length_top, width_left, width_right, length_bottom );
    break;
  case PLA_MATRIX:
    PLA_Obj_view_shift_matrix( 
	  obj, length_top, width_left, width_right, length_bottom );
    break;
  case PLA_MSCALAR:
    PLA_Obj_view_shift_mscalar( 
	  obj, length_top, width_left, width_right, length_bottom );
    break;
  case PLA_PMVECTOR:
    PLA_Obj_view_shift_pmvector( 
	  obj, length_top, width_left, width_right, length_bottom );
    break; 
  }

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    return_value = PLA_Obj_view_shift_exit( 
	  obj, length_top, width_left, width_right, length_bottom );

  return return_value;
}

int PLA_Obj_view_shift_mv (PLA_Obj      obj,           int      length_top,
                          int          width_left,    int      width_right,
                          int          length_bottom)

/*----------------------------------------------------------------------------

Purpose : Shift the boundaries of the view of the linear algebra object



IN/OUT       obj           object to be grown [shifted]
IN           length_top    length to shift top boundary
IN           width_left    width to shift left boundary
IN           width_right   width to shift right boundary
IN           length_bottom length to shift bottom boundary

Note : Positive is down/right negative is up/left (as if you were in quadrant IV
in the cartesian plane.

----------------------------------------------------------------------------*/
{
  PLA_Obj 
    view;

  PLA_Template
    templ;

  int
    diff_local_length, diff_local_width, local_ldim, typesize,
    zero_or_one,
    value = PLA_SUCCESS;

  void 
    *local_buffer;
  
  MPI_Datatype
    datatype;

  view = obj;

  if ( length_top == 0 ) diff_local_length = 0;
  else if ( length_top > 0 ) 
    PLA_Temp_compute_local_length( PLA_COMM_ALL, view->base_obj->templ,
				    view->global_align_row, length_top,
				    &diff_local_length );
  else {
    PLA_Temp_compute_local_length( PLA_COMM_ALL, view->base_obj->templ,
				    view->global_align_row+length_top, -length_top,
				    &diff_local_length );
    diff_local_length = - diff_local_length;
  }

  diff_local_width = width_left;

  PLA_Obj_datatype( obj, &datatype );
  MPI_Type_size( datatype, &typesize );

  PLA_Obj_template( obj, &templ );
    
  view->global_length     += length_bottom - length_top;
  view->global_width      += width_right - width_left;
  view->global_align_row  += length_top;

  PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ALL,                   /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       view->global_align_row,        /* Alignment to template wrt direction */
       view->global_length,           /* Length in direction */ 
       &view->local_length,           /* Length assigned to this node */
       &view->size_top,               /* Size of first block */
       &view->owner_top,	      /* Owner of first block */
       &view->size_bottom,            /* Size of last block */
       &view->owner_bottom );         /* Owner of last block */

  view->size_left         = view->global_width;
  view->size_right        = view->global_width;
    
  view->local_width       = view->global_width;

  PLA_Obj_local_ldim( obj, &local_ldim );
  PLA_Obj_local_buffer( obj, &local_buffer );
  view->local_buffer = 
    (char *) local_buffer +
    ( diff_local_width * local_ldim + diff_local_length ) * typesize;

  return value;
}

int PLA_Obj_view_shift_matrix
                         (PLA_Obj      obj,           int      length_top,
                          int          width_left,    int      width_right,
                          int          length_bottom)

/*----------------------------------------------------------------------------

Purpose : Shift the boundaries of the view of the linear algebra object



IN/OUT       obj           object to be grown [shifted]
IN           length_top    length to shift top boundary
IN           width_left    width to shift left boundary
IN           width_right   width to shift right boundary
IN           length_bottom length to shift bottom boundary

Note : Positive is down/right negative is up/left (as if you were in quadrant IV
in the cartesian plane.

----------------------------------------------------------------------------*/
{
  PLA_Obj 
    view;

  PLA_Template
    templ;

  int
    zero_or_one,
    value = PLA_SUCCESS,
    diff_local_length, diff_local_width, local_ldim, typesize;

  MPI_Datatype
    datatype;

  void 
    *local_buffer;

  PLA_Obj_template( obj, &templ );
  view     = obj;

  if ( length_top == 0 ) diff_local_length = 0;
  else if ( length_top > 0 ) 
    PLA_Temp_compute_local_length( PLA_COMM_COL, view->base_obj->templ,
				    view->global_align_row, length_top,
				    &diff_local_length );
  else {
    PLA_Temp_compute_local_length( PLA_COMM_COL, view->base_obj->templ,
				    view->global_align_row+length_top, -length_top,
				    &diff_local_length );
    diff_local_length = - diff_local_length;
  }

  if ( width_left == 0 ) diff_local_width = 0;
  else if ( width_left > 0 ) 
    PLA_Temp_compute_local_length( PLA_COMM_ROW, view->base_obj->templ,
				    view->global_align_col, width_left,
				    &diff_local_width );
  else {
    PLA_Temp_compute_local_length( PLA_COMM_ROW, view->base_obj->templ,
				    view->global_align_col+width_left, -width_left,
				    &diff_local_width );
    diff_local_width = - diff_local_width;
  }

  PLA_Obj_datatype( obj, &datatype );
  MPI_Type_size( datatype, &typesize );

  view->global_length     += length_bottom - length_top;
  view->global_width      += width_right - width_left;
  view->global_align_row  += length_top;
  view->global_align_col  += width_left;

  PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_COL,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       view->global_align_row,        /* Alignment to template wrt direction */
       view->global_length,                 /* Length in direction */ 
       &view->local_length,           /* Length assigned to this node */
       &view->size_top,               /* Size of first block */
       &view->owner_top,	      /* Owner of first block */
       &view->size_bottom,            /* Size of last block */
       &view->owner_bottom );         /* Owner of last block */

  PLA_Temp_compute_owners_and_sizes(
       PLA_COMM_ROW,                 /* Direction in which to compute */
       templ,                         /* Template for distribution */   
       view->global_align_col,        /* Alignment to template wrt direction */
       view->global_width,                  /* Length in direction */ 
       &view->local_width,            /* Length assigned to this node */
       &view->size_left,              /* Size of first block */
       &view->owner_left,             /* Owner of first block */
       &view->size_right,             /* Size of last block */
       &view->owner_right );          /* Owner of last block */

  PLA_Obj_local_ldim( obj, &local_ldim );
  PLA_Obj_local_buffer( obj, &local_buffer );
  view->local_buffer = 
    (char *) local_buffer +
    ( diff_local_width * local_ldim + diff_local_length ) * typesize;
 
  return value;
}

int PLA_Obj_view_shift_mscalar
                         (PLA_Obj      obj,           int      length_top,
                          int          width_left,    int      width_right,
                          int          length_bottom)

/*----------------------------------------------------------------------------

Purpose : Shift the boundaries of the view of the linear algebra object



IN/OUT       obj           object to be grown [shifted]
IN           length_top    length to shift top boundary
IN           width_left    width to shift left boundary
IN           width_right   width to shift right boundary
IN           length_bottom length to shift bottom boundary

Note : Positive is down/right negative is up/left (as if you were in quadrant IV
in the cartesian plane.

----------------------------------------------------------------------------*/
{
  PLA_Obj 
    view;

  PLA_Template
    templ;
  int
    diff_local_length, diff_local_width, local_ldim, typesize,
    myrow, mycol,
    value = PLA_SUCCESS;

  void 
    *local_buffer;

  MPI_Datatype
    datatype;
  
  view     = obj;
    
  view->global_length     += length_bottom - length_top;
  view->global_width      += width_right - width_left;

  PLA_Obj_template( obj, &templ );
  PLA_Temp_comm_col_rank( templ, &myrow );
  PLA_Temp_comm_row_rank( templ, &mycol )
;
  /* Change this later! */
  if ( ( view->owner_row == PLA_ALL_ROWS || view->owner_row == myrow ) &&
       ( view->owner_col == PLA_ALL_COLS || view->owner_col == mycol ) ) {
    view->local_length += length_bottom - length_top;
    view->local_width  += width_right - width_left;
  }

  view->size_top  = view->size_bottom = view->global_length;
  view->size_left = view->size_right  = view->global_width;

  PLA_Obj_datatype( obj, &datatype );
  MPI_Type_size( datatype, &typesize );

  diff_local_length = length_top;
  diff_local_width  = width_left;

  PLA_Obj_local_ldim( obj, &local_ldim );
  PLA_Obj_local_buffer( obj, &local_buffer );
  view->local_buffer = 
    (char *) local_buffer +
    ( diff_local_width * local_ldim + diff_local_length ) * typesize;

  return value;
}

int PLA_Obj_view_shift_pmvector
                         (PLA_Obj      obj,           int      length_top,
                          int          width_left,    int      width_right,
                          int          length_bottom)

/*----------------------------------------------------------------------------

Purpose : Shift the boundaries of the view of the linear algebra object



IN/OUT       obj           object to be grown [shifted]
IN           length_top    length to shift top boundary
IN           width_left    width to shift left boundary
IN           width_right   width to shift right boundary
IN           length_bottom length to shift bottom boundary

Note : Positive is down/right negative is up/left (as if you were in quadrant IV
in the cartesian plane.

----------------------------------------------------------------------------*/
{
  PLA_Obj_base 
    base_obj;

  PLA_Obj 
    view;

  PLA_Template
    templ;

  int
    proj_onto, diff_local_length, diff_local_width, typesize, local_ldim,
    owner,
    project_onto,
    value = PLA_SUCCESS,
    myrow, mycol;

  void
    *local_buffer;

  MPI_Datatype
    datatype;

  view     = obj;

  PLA_Obj_project_onto( obj, &proj_onto );

  if ( proj_onto == PLA_PROJ_ONTO_COL ){
    if ( length_top == 0 ) diff_local_length = 0;
    else if ( length_top > 0 ) 
      PLA_Temp_compute_local_length( PLA_COMM_COL, view->base_obj->templ,
				      view->global_align_row, length_top,
				      &diff_local_length );
    else {
      PLA_Temp_compute_local_length( PLA_COMM_COL, view->base_obj->templ,
				      view->global_align_row+length_top, -length_top,
				      &diff_local_length );
      diff_local_length = - diff_local_length;
    }
    diff_local_width = width_left;
  }
  else{
    if ( width_left == 0 ) diff_local_width = 0;
    else if ( width_left > 0 ) 
      PLA_Temp_compute_local_length( PLA_COMM_ROW, view->base_obj->templ,
				      view->global_align_col, width_left,
				      &diff_local_width );
    else {
      PLA_Temp_compute_local_length( PLA_COMM_ROW, view->base_obj->templ,
				      view->global_align_col+width_left, -width_left,
				      &diff_local_width );
      diff_local_width = - diff_local_width;
    }

    diff_local_length = length_top;
  }

  PLA_Obj_datatype( obj, &datatype );
  MPI_Type_size( datatype, &typesize );

  PLA_Obj_template( obj, &templ );
  view->global_length     += length_bottom - length_top;
  view->global_width      += width_right - width_left;
  PLA_Obj_project_onto( obj, &project_onto );
  if ( project_onto == PLA_PROJ_ONTO_COL )
    view->global_align_row  += length_top;
  else
    view->global_align_col  += width_left;

  if ( project_onto == PLA_PROJ_ONTO_COL ){
    PLA_Temp_compute_owners_and_sizes(
          PLA_COMM_COL,              /* Direction in which to compute */
          templ,                      /* Template for distribution */   
          view->global_align_row,     /* Alignment to template wrt direction */
          view->global_length,         /* Length in direction */ 
          &view->local_length,        /* Length assigned to this node */
          &view->size_top,            /* Size of first block */
          &view->owner_top,	      /* Owner of first block */
          &view->size_bottom,         /* Size of last block */
          &view->owner_bottom );      /* Owner of last block */

    owner = obj->owner_col;
    PLA_Temp_comm_row_rank( templ, &mycol );
    if ( mycol == owner || PLA_ALL_COLS == owner )
      view->local_width      = view->global_width;
    else 
      view->local_width = 0;

    view->size_left         = view->global_width;
    view->size_right        = view->global_width;
  }
  else {
    PLA_Temp_compute_owners_and_sizes(
          PLA_COMM_ROW,              /* Direction in which to compute */
          templ,                      /* Template for distribution */   
          view->global_align_col,     /* Alignment to template wrt direction */
          view->global_width,          /* Length in direction */ 
          &view->local_width,         /* Length assigned to this node */
          &view->size_left,           /* Size of first block */
          &view->owner_left,	      /* Owner of first block */
          &view->size_right,          /* Size of last block */
          &view->owner_right );       /* Owner of last block */

    PLA_Temp_comm_col_rank( templ, &myrow );
    owner = obj->owner_row;
    if ( myrow == owner || PLA_ALL_ROWS == owner )
      view->local_length   = view->global_length;
    else 
      view->local_length = 0;

    view->size_top         = view->global_length;
    view->size_bottom      = view->global_length;
  }

  PLA_Obj_local_ldim( obj, &local_ldim );
  PLA_Obj_local_buffer( obj, &local_buffer );
  view->local_buffer = 
    (char *) local_buffer +
    ( diff_local_width * local_ldim + diff_local_length ) * typesize;

  return value;
}

/***************************************************************************/

int PLA_Obj_split_size   (PLA_Obj      obj,      int      side,
                           int          *size,    int      *owner)

/*----------------------------------------------------------------------------

Purpose : Compute size of split to next block boundary.

IN       obj       object to be split
IN       side      side of split
OUT      size      size to next template subblock split
OUT      owner     index of row or column of nodes that owns the split block

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_split_size_enter( obj, side, size, owner );

  if ( side == PLA_SIDE_TOP ){
    PLA_OBJ_SIZE_TOP     ( obj, size );
    PLA_OBJ_OWNER_TOP    ( obj, owner );
  }
  else if ( side == PLA_SIDE_BOTTOM ){
    PLA_OBJ_SIZE_BOTTOM  ( obj, size );
    PLA_OBJ_OWNER_BOTTOM ( obj, owner );
  }
  else if ( side == PLA_SIDE_LEFT ){
    PLA_OBJ_SIZE_LEFT    ( obj, size );
    PLA_OBJ_OWNER_LEFT   ( obj, owner );
  }
  else if ( side == PLA_SIDE_RIGHT ){
    PLA_OBJ_SIZE_RIGHT   ( obj, size );
    PLA_OBJ_OWNER_RIGHT  ( obj, owner );
  }
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_split_size_exit( obj, side, size, owner );

  return value;
}

/***************************************************************************/

int PLA_Vector_create_conf_to   (PLA_Obj      obj,        PLA_Obj     *new_vector)

/*----------------------------------------------------------------------------

Purpose : Create vector conformal to given object.


IN      obj          original object
OUT     new_vector   created object

-------------------------------------------------------------------------*/
{
  PLA_Abort( "PLA_Vector_create_conf_to: not supported", __LINE__, __FILE__ );
}

/***************************************************************************/

int PLA_Mvector_create_conf_to   (PLA_Obj      obj,        int      global_width,
                                  PLA_Obj     *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create mvector conformal to given object.


IN      obj                original object
IN      global_width       number of vectors in mvector
OUT     new_obj            created object

-------------------------------------------------------------------------*/
{
  PLA_Obj_base 
    base_obj;

  PLA_Obj 
    view;

  PLA_Template
    templ;

  int
    zero_or_one,
    value = PLA_SUCCESS,
    objtype,
    typesize,
    i, 
    global_length, 
    global_align_row, 
    project_onto;

  void 
    *buffer = NULL;
  
  char 
    *tempp_new, *tempp_old;
  
  MPI_Datatype 
    datatype;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Mvector_create_conf_to_enter( obj, global_width, new_obj );

  PLA_OBJ_OBJTYPE( obj, &objtype );
  PLA_OBJ_DATATYPE( obj, &datatype );
  PLA_OBJ_PROJECT_ONTO( obj, &project_onto );
  if ( project_onto == PLA_PROJ_ONTO_COL ){
    PLA_OBJ_GLOBAL_LENGTH( obj, &global_length );
    PLA_OBJ_GLOBAL_ALIGN_ROW( obj, &global_align_row );
  }
  else {
    PLA_OBJ_GLOBAL_WIDTH( obj, &global_length );
    PLA_OBJ_GLOBAL_ALIGN_COL( obj, &global_align_row );
  }
  PLA_OBJ_TEMPLATE( obj, &templ );
  PLA_Temp_zero_or_one( templ, &zero_or_one );
  
  switch( objtype ){
  case PLA_MATRIX:
  case PLA_PMVECTOR: {
    value = PLA_Mvector_create( datatype, global_length, global_width, templ, 
				  global_align_row, new_obj );
    break;
  }
  case PLA_MVECTOR: {
    if ( *new_obj != NULL )
      PLA_Obj_free( new_obj );
    
    base_obj = pla_get_base_obj( );
    view     = pla_get_view( );
    
    /* copy contents of current object */
    tempp_new = ( char * ) view;
    tempp_old = ( char * ) obj;
    
    for ( i=0; i<sizeof( struct PLA_Obj_view_struct ); i++ )
      *tempp_new++ = *tempp_old++;
    
    base_obj->datatype         = datatype;
    
    base_obj->global_length    = global_length;
    base_obj->global_width     = global_width;
    
    base_obj->global_align_row = global_align_row;
    base_obj->global_align_col = zero_or_one;
    
    base_obj->local_buffer     = NULL;
    base_obj->user_buffer      = PLA_UNDEFINED;
    base_obj->number_of_views  = 1;
    base_obj->local_buffer_size = 0;
    base_obj->ooc              = FALSE;
    
    base_obj->templ            = templ;
    
#if DEBUG==1
    view->pla_obj           = NULL;
#endif

    view->fixed             = FALSE; 
    view->global_width      = global_width;
    
    view->global_align_col  = zero_or_one;
    
    view->size_left         = global_width;
    view->size_right        = global_width;
    
    view->local_width       = global_width;
    view->local_buffer      = NULL;
    base_obj->local_ldim    = view->local_length;
    
    view->base_obj          = base_obj;

    view->mode   = PLA_MODE_CLOSED;

    MPI_Type_size( datatype, &typesize );
    if ( view->local_length * global_width != 0 ){
      buffer = PLA_malloc( (size_t) typesize * view->local_length * global_width );
      
      if ( buffer == NULL ) 
	PLA_Abort( "PLA_Mvector_create_conf_to: malloc failed\n", __LINE__, __FILE__ );
    }
    
    PLA_Obj_attach_buffer( view, buffer, view->local_length, FALSE );
    
    *new_obj = view;

    PLA_Obj_set_to_zero( *new_obj );
    break;
  }
  default: {
    value--;
  }
  }
  
  if ( PLA_ERROR_CHECKING )    /* Perform exit error checking */
    value = PLA_Mvector_create_conf_to_exit( obj, global_width, new_obj );

  return value;
}

/***************************************************************************/

int PLA_Pvector_create_conf_to   (PLA_Obj      obj,        int         project_onto,
                                  int          owner,      PLA_Obj     *new_pvector)

/*----------------------------------------------------------------------------

Purpose : Create projected vector conformal to given object.


IN      obj                original object
IN      project_onto       mesh direction onto which to project
IN      owner              index of row or column of nodes which contains owning nodes
OUT     new_pvector        created object

-------------------------------------------------------------------------*/
{
  PLA_Abort( "PLA_Pvector_create_conf_to: not supported", __LINE__, __FILE__ );
}

/***************************************************************************/

int PLA_Pmvector_create_conf_to   (PLA_Obj      obj,        int         project_onto,
                                   int          owner,      int         num_vectors,
                                   PLA_Obj     *new_pmvector)

/*----------------------------------------------------------------------------

Purpose : Create projected mvector conformal to given object.


IN      obj                original object
IN      project_onto       mesh direction onto which to project
IN      owner              index of row or column of nodes which contains owning nodes
IN      num_vectors        number of vectors in multivector
OUT     new_pmvector       created object

-------------------------------------------------------------------------*/
{
  int 
    value,
    objtype,
    length, width,
    proj_onto_from,
    align, align_row, align_col;
  MPI_Datatype
    datatype;
  PLA_Template
    templ;

  PLA_Obj_objtype( obj, &objtype );
  PLA_Obj_datatype( obj, &datatype );
  PLA_Obj_template( obj, &templ );

  switch( objtype ){
  case PLA_MATRIX:
    switch( project_onto ){
    case PLA_PROJ_ONTO_COL:
      PLA_Obj_global_length( obj, &length );
      PLA_Obj_global_align_row( obj, &align_row );
      value = PLA_Pmvector_create( datatype, project_onto, owner, length,
				    num_vectors, templ, align_row, 
				    new_pmvector );
      break;
    case PLA_PROJ_ONTO_ROW:
      PLA_Obj_global_width( obj, &width );
      PLA_Obj_global_align_col( obj, &align_col );
      value = PLA_Pmvector_create( datatype, project_onto, owner, num_vectors,
				    width, templ, align_col, new_pmvector );
      break;
    }
    break;
  case PLA_MVECTOR:
    PLA_Obj_global_length( obj, &length );
    PLA_Obj_global_align_row( obj, &align );
    switch( project_onto ){
    case PLA_PROJ_ONTO_COL:
      value = PLA_Pmvector_create( datatype, project_onto, owner, length,
				    num_vectors, templ, align, 
				    new_pmvector );
      break;
    case PLA_PROJ_ONTO_ROW:
      value = PLA_Pmvector_create( datatype, project_onto, owner, num_vectors,
				    length, templ, align, new_pmvector );
      break;
    }
    break;
  case PLA_PMVECTOR:
    PLA_Obj_project_onto( obj, &proj_onto_from );
    switch( proj_onto_from ){
    case PLA_PROJ_ONTO_COL:
      PLA_Obj_global_length( obj, &length );
      PLA_Obj_global_align_row( obj, &align );
      break;
    case PLA_PROJ_ONTO_ROW:
      PLA_Obj_global_width( obj, &length );
      PLA_Obj_global_align_col( obj, &align );
      break;
    }

    switch( project_onto ){
    case PLA_PROJ_ONTO_COL:
      value = PLA_Pmvector_create( datatype, project_onto, owner, length,
				    num_vectors, templ, align, 
				    new_pmvector );
      break;
    case PLA_PROJ_ONTO_ROW:
      value = PLA_Pmvector_create( datatype, project_onto, owner, num_vectors,
				    length, templ, align, new_pmvector );
      break;
    }
    break;
  default: 
    PLA_Abort( "PLA_Pmvector_create_conf_to: objtype not yet supported \n", __LINE__, __FILE__ );
  }
  
  return value;
}

 

/***************************************************************************/

int PLA_Matrix_create_conf_to   (PLA_Obj      matrix,        PLA_Obj     *new_matrix)

/*----------------------------------------------------------------------------

Purpose : Create matrix conformal to given object (must be a matrix).



IN      matrix       original object
OUT     new_matrix   created object

-------------------------------------------------------------------------*/
{
  int 
    objtype,
    length, width,
    align_row, align_col;
  MPI_Datatype
    datatype;
  PLA_Template 
    templ = NULL;
  
  /* OPTIMIZE LATER!!! */
  PLA_Obj_objtype( matrix, &objtype );
  if ( objtype != PLA_MATRIX )
    PLA_Abort( "illegal objtype in PLA_Matrix_create_conf_to", __LINE__, __FILE__ );

  PLA_Obj_datatype( matrix, &datatype );
  PLA_Obj_global_length( matrix, &length );
  PLA_Obj_global_width( matrix, &width );
  PLA_Obj_global_align_row( matrix, &align_row );
  PLA_Obj_global_align_col( matrix, &align_col );
  PLA_Obj_template( matrix, &templ );
  
  return PLA_Matrix_create( datatype, length, width, templ, align_row, 
		      align_col, new_matrix );
}

/***************************************************************************/

int PLA_Mscalar_create_conf_to   (PLA_Obj      obj,        int         owner_row,
                                  int          owner_col,  PLA_Obj     *new_mscalar)

/*----------------------------------------------------------------------------

Purpose : Create multiscalar conformal to given object.


IN      obj           original object
IN      owner_row     index of row of nodes that contains owning node(s)
IN      owner_col     index of column of nodes that contains owning node(s)
OUT     new_mscalar   created object (multiscalar)

----------------------------------------------------------------------------*/
{
  int 
    value,
    objtype,
    length, width;
  MPI_Datatype
    datatype;
  PLA_Template 
    templ = NULL;
  
  PLA_Obj_datatype( obj, &datatype );
  PLA_Obj_global_length( obj, &length );
  PLA_Obj_global_width( obj, &width );
  PLA_Obj_template( obj, &templ );

  if ( owner_row == PLA_INHERIT || owner_col == PLA_INHERIT )
    PLA_Abort("PLA_Mscalar_create_conf_to: PLA_INHERIT not yet supported", __LINE__, __FILE__ );

  value = PLA_Mscalar_create( datatype, owner_row, owner_col, 
			      length, width, templ, new_mscalar );

  return value;
}

/***************************************************************************/

int PLA_Obj_objtype_cast   (PLA_Obj      obj,         int      obj_type)

/*----------------------------------------------------------------------------

Purpose : Cast the linear algebra object to have given object type.


IN/OUT        obj            original object
IN            objtype        new object type

-------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    objtype_old;
  PLA_Template
    templ;
  int 
    nprows;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_objtype_cast_enter( obj, obj_type );

  if ( value == 0 ){
    PLA_Obj_objtype( obj, &objtype_old );
    switch( objtype_old ){
    case PLA_MATRIX:
      switch( obj_type ){
      case PLA_MSCALAR:
	obj->owner_row = obj->owner_top;
	obj->owner_col = obj->owner_left;
	obj->objtype   = obj_type;
	break;
      default:
	PLA_Warning( "objtype of cast not yet supported" );
	break;
      }
      break;
    case PLA_MVECTOR:
      switch( obj_type ){
      case PLA_MSCALAR:
	PLA_Obj_template( obj, &templ );
	PLA_Temp_comm_col_size( templ, &nprows );
	obj->owner_row = obj->owner_top%nprows;
	obj->owner_col = obj->owner_top/nprows;
	obj->objtype   = obj_type;
	break;
      default:
	PLA_Warning( "objtype of cast not yet supported" );
	break;
      }
      break;
    default:
      PLA_Warning( "objtype not yet supported" );
      break;
    }
  }


  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_objtype_cast_exit( obj, obj_type );

  return value;
}

/***************************************************************************/

int PLA_Obj_set_orientation   (PLA_Obj     obj,     int     project_onto)

/*----------------------------------------------------------------------------

Purpose : Set the orientation of a matrix panel.


IN/OUT      obj              original object
IN          project_onto     annotation indicating orientation

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_set_orientation_enter( obj, project_onto );

  obj->proj_onto = project_onto;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_set_orientation_exit( obj, project_onto );

  return value;
}

/***************************************************************************/

int PLA_Obj_get_orientation   (PLA_Obj     obj,     int     *project_onto)


/*----------------------------------------------------------------------------

Purpose : Return the orientation of a matrix panel.


IN/OUT      obj              original object
IN          project_onto     annotation indicating orientation

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  
  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_get_orientation_enter( obj, project_onto );

  PLA_OBJ_PROJECT_ONTO( obj, project_onto );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Obj_get_orientation_exit( obj, project_onto );

  return value;
}

/***************************************************************************/

int PLA_Create_constants_conf_to( PLA_Obj obj, 
                                  PLA_Obj *minus_one, 
                                  PLA_Obj *zero,
                                  PLA_Obj *one )

/*----------------------------------------------------------------------------

Purpose : Create duplicated constants -1, 0, and/or 1 of the appropriate type
          and template, as dictated by obj.

IN          obj              original object
OUT         minus_one        object containing duplicated multiscalar -1
OUT         zero             object containing duplicated multiscalar 0
OUT         one              object containing duplicated multiscalar 1

if any of minus_one, zero, or one come in as NULL, no action is taken
for that parameter.
----------------------------------------------------------------------------*/
{
  PLA_Template templ;
  MPI_Datatype datatype;

  PLA_Obj_datatype( obj, &datatype );
  PLA_Obj_template( obj, &templ );

  if ( one != NULL ) {
    PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
                        1, 1, templ, one );
    PLA_Obj_set_to_one( *one );
  }

  if ( zero != NULL ) {
    PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
                        1, 1, templ, zero );
    PLA_Obj_set_to_zero( *zero );
  }

  if ( minus_one != NULL ) {
    PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS,
		       1, 1, templ, minus_one );
    PLA_Obj_set_to_minus_one( *minus_one );
  }

  return( PLA_SUCCESS );
}

