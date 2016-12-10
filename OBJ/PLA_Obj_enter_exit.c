/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/***************************************************************************/

int PLA_Vector_create_enter   (
           MPI_Datatype   datatype,   int    global_length,
	   PLA_Template   templ,      int    global_align,
	   PLA_Obj       *new_obj)

/*----------------------------------------------------------------------------

Purpose : Create distributed vector.


IN     datatype          datatype of object
IN     global_length     global length of vector
IN     templ             template for vector and matrix distribution
IN     global_align      alignment to template
OUT    new_obj           object describing created vector

----------------------------------------------------------------------------*/
{
}


int PLA_Vector_create_exit   (
           MPI_Datatype   datatype,   int    global_length,
	   PLA_Template   templ,      int    global_align,
	   PLA_Obj       *new_obj)

/*----------------------------------------------------------------------------

Purpose : Create distributed vector.


IN     datatype          datatype of object
IN     global_length     global length of vector
IN     templ             template for vector and matrix distribution
IN     global_align      alignment to template
OUT    new_obj           object describing created vector

----------------------------------------------------------------------------*/
{
}

/***************************************************************************/

int PLA_Mvector_create_without_buffer_enter (
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
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Mvector_create_wout_buffer";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( !PLA_Valid_datatype( datatype ) ){
      PLA_Warning( "Invalid datatype" );
      value--;
    }
    if ( global_length < 0 ){
      PLA_Warning( "Negative global length" );
      value--;
    }
    if ( global_width < 0 ){
      PLA_Warning( "Negative global width" );
      value--;
    }
    if ( !PLA_Valid_template( templ ) ){
      PLA_Warning( "Invalid template" );
      value--;
    }
    if ( global_align < 0 && global_align != PLA_ALIGN_FIRST ){
      PLA_Warning( "Negative global align" );
      value--;
    }
    if ( *new_obj != NULL && !PLA_Valid_object( *new_obj ) ) {
      PLA_Warning( "Invalid output object parameter" );
      value--;
    }
  }

  return value;
}


int PLA_Mvector_create_without_buffer_exit (
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
  int
    value = PLA_SUCCESS;

  char 
    routine_name[ 35];

#if DEBUG==1  
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Mvector_create( datatype, global_length, global_width,
                        templ, 
		        ( global_align == PLA_ALIGN_FIRST ? 
			 PLA_ALIGN_FIRST : global_align ), 
		        &(*new_obj)->pla_obj );

    if ( !PLA_Obj_consistent_with_R12( *new_obj ) ){
      PLA_Warning( "Objects inconsistent" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Mvector_create_enter   (
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
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Mvector_create_exit   (
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

  return value;
}

/***************************************************************************/

int PLA_Matrix_create_without_buffer_enter   (
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
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Matrix_create_wout_buffer";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( !PLA_Valid_datatype( datatype ) ){
      PLA_Warning( "Invalid datatype" );
      value--;
    }
    if ( global_length < 0 ){
      PLA_Warning( "Negative global length" );
      value--;
    }
    if ( global_width < 0 ){
      PLA_Warning( "Negative global width" );
      value--;
    }
    if ( !PLA_Valid_template( templ ) ){
      PLA_Warning( "Invalid template" );
      value--;
    }
    if ( global_align_row < 0 && global_align_row != PLA_ALIGN_FIRST  ){
      PLA_Warning( "Negative global row align" );
      value--;
    }

    if ( global_align_col < 0 && global_align_col != PLA_ALIGN_FIRST  ){
      PLA_Warning( "Negative global col align" );
      value--;
    }

    if ( *new_obj != NULL && !PLA_Valid_object( *new_obj ) ) {
      PLA_Warning( "Invalid output object parameter" );
      value--;
    }
  }

  return value;
}


int PLA_Matrix_create_without_buffer_exit   (
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
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Matrix_create( datatype, global_length, global_width,
          templ, 
          ( global_align_row == PLA_ALIGN_FIRST ? 
	    PLA_ALIGN_FIRST : global_align_row ), 
          ( global_align_col == PLA_ALIGN_FIRST ? 
	    PLA_ALIGN_FIRST : global_align_col ), 
	  &(*new_obj)->pla_obj );

    if ( !PLA_Obj_consistent_with_R12( *new_obj ) ){
      PLA_Warning( "Objects inconsistent" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Matrix_create_enter   (
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
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Matrix_create_exit   (
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
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Mscalar_create_without_buffer_enter  (
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
  int 
    value = PLA_SUCCESS,
    nprows, npcols;
  char 
    routine_name[ 35 ] = "PLA_Mscalar_create_wt_buffer";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( !PLA_Valid_datatype( datatype ) ){
      PLA_Warning( "Invalid datatype" );
      value--;
    }

    PLA_Temp_comm_col_size( templ, &nprows );
    PLA_Temp_comm_row_size( templ, &npcols );

    if ( ( owner_row < 0 || owner_row >= nprows ) &&
	   owner_row != PLA_ALL_ROWS && owner_row != PLA_UNDEFINED ) {
      PLA_Warning( "Illegal owner_row" );
      value--;
    }
    if ( ( owner_col < 0 || owner_col >= npcols ) &&
	   owner_col != PLA_ALL_COLS ) {
      PLA_Warning( "Illegal owner_col" );
      value--;
    }
    if ( length < 0 ) {
      PLA_Warning( "Negative global length" );
      value--;
    }
    if ( width < 0 ){
      PLA_Warning( "Negative global width" );
      value--;
    }
    if ( !PLA_Valid_template( templ ) ){
      PLA_Warning( "Invalid template" );
      value--;
    }
    if ( *new_obj != NULL && !PLA_Valid_object( *new_obj ) ) {
      PLA_Warning( "Invalid output object parameter" );
      value--;
    }
  }

  return value;
}


int PLA_Mscalar_create_without_buffer_exit  (
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
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Mscalar_create( 
          datatype, 
	  ( owner_row == PLA_ALL_ROWS ? PLA_ALL_ROWS : owner_row ),
	  ( owner_col == PLA_ALL_COLS ? PLA_ALL_COLS : owner_col ),
          length, width, templ, 
          &(*new_obj)->pla_obj );

    if ( !PLA_Obj_consistent_with_R12( *new_obj ) ){
      PLA_Warning( "Objects inconsistent" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Mscalar_create_enter  (
       MPI_Datatype   datatype,          int            owner_row,
       int            owner_col,         int            length,
       int            width,             PLA_Template   templ,
       PLA_Obj       *new_mscalar)

/*----------------------------------------------------------------------------

Purpose : Create distributed mutltiscalar.


IN     datatype          datatype of object
IN     owner_row         index of row of nodes with contains owning node(s)
IN     owner_col         index of column of nodes with contains owning node(s)
IN     length            length of mscalar
IN     width             width of mscalar
IN     templ             template for vector and matrix distribution
OUT    new_mscalar       object describing created multiscalar

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Mscalar_create_exit  (
       MPI_Datatype   datatype,          int            owner_row,
       int            owner_col,         int            length,
       int            width,             PLA_Template   templ,
       PLA_Obj       *new_mscalar)

/*----------------------------------------------------------------------------

Purpose : Create distributed mutltiscalar.


IN     datatype          datatype of object
IN     owner_row         index of row of nodes with contains owning node(s)
IN     owner_col         index of column of nodes with contains owning node(s)
IN     length            length of mscalar
IN     width             width of mscalar
IN     templ             template for vector and matrix distribution
OUT    new_mscalar       object describing created multiscalar

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Pvector_create_enter  (
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
  int value = PLA_SUCCESS;

  return value;
}



int PLA_Pvector_create_exit  (
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
  int value = PLA_SUCCESS;

  return value;
}


/***************************************************************************/

int PLA_Pmvector_create_without_buffer_enter  (
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
IN     templ                 template for vector and matrix distribution
IN     global_align          alignment to template
OUT    new_obj               object describing created projected multi-vector

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    nprows, npcols;
  char 
    routine_name[ 35 ] = "PLA_Pmvector_create_wout_buffer";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( !PLA_Valid_datatype( datatype ) ){
      PLA_Warning( "Invalid datatype" );
      value--;
    }
    if ( project_onto != PLA_PROJ_ONTO_ROW && 
         project_onto != PLA_PROJ_ONTO_COL ){
      PLA_Warning( "Illegal value for project_onto" );
      value--;
    }
    if ( project_onto == PLA_PROJ_ONTO_ROW ){
      PLA_Temp_comm_col_size( templ, &nprows );
      if ( owner != PLA_ALL_ROWS && ( owner < 0 || owner >= nprows ) ){
	PLA_Warning( "Illegal owner" );
	value--;
      }
    }
    else if ( project_onto == PLA_PROJ_ONTO_COL ){
      PLA_Temp_comm_row_size( templ, &npcols );
      if ( owner != PLA_ALL_COLS && ( owner < 0 || owner >= npcols ) ){
	PLA_Warning( "Illegal owner" );
	value--;
      }
    }
    if ( global_proj_width < 0 ){
      PLA_Warning( "Negative global width" );
      value--;
    }
    if ( !PLA_Valid_template( templ ) ){
      PLA_Warning( "Invalid template" );
      value--;
    }
    if ( global_align < 0 && global_align != PLA_ALIGN_FIRST ){
      PLA_Warning( "Negative global align" );
      value--;
    }
    if ( *new_obj != NULL && !PLA_Valid_object( *new_obj ) ) {
      PLA_Warning( "Invalid output object parameter" );
      value--;
    }
  }

  return value;
}


int PLA_Pmvector_create_without_buffer_exit  (
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
IN     templ                 template for vector and matrix distribution
IN     global_align          alignment to template
OUT    new_obj               object describing created projected multi-vector

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Pmvector_create( 
          datatype, 
          ( project_onto == PLA_PROJ_ONTO_ROW ? 
	        PLA_PROJ_ONTO_ROW : PLA_PROJ_ONTO_COL), 
	  ( owner == PLA_ALL_ROWS ? PLA_ALL_ROWS : 
               ( owner == PLA_ALL_COLS ? PLA_ALL_COLS : owner ) ),
	  global_proj_length, 
          global_proj_width, templ, 
	  (global_align == PLA_ALIGN_FIRST ? PLA_ALIGN_FIRST : global_align),
          &(*new_obj)->pla_obj );

    if ( !PLA_Obj_consistent_with_R12( *new_obj ) ){
      PLA_Warning( "Objects inconsistent" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Pmvector_create_enter  (
       MPI_Datatype   datatype,           int            project_onto,
       int            owner,              int            global_proj_length,
       int            global_proj_width,  PLA_Template   templ,
       int            global_align,       PLA_Obj       *new_proj_mvector)

/*----------------------------------------------------------------------------

Purpose : Create projected multi-vector.

IN     datatype              datatype of object
IN     project_onto          direction onto which to project
IN     owner                 index of row or column of nodes wich containts projected vector
IN     global_proj_length    length of (projected) vector
IN     global_proj_width     width of (projected) vector
IN     global_align          alignment to template
IN     templ                 template for vector and matrix distribution
OUT    new_proj_mvector      object describing created projected multi-vector

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Pmvector_create_exit  (
       MPI_Datatype   datatype,           int            project_onto,
       int            owner,              int            global_proj_length,
       int            global_proj_width,  PLA_Template   templ,
       int            global_align,       PLA_Obj       *new_proj_mvector)

/*----------------------------------------------------------------------------

Purpose : Create projected multi-vector.

IN     datatype              datatype of object
IN     project_onto          direction onto which to project
IN     owner                 index of row or column of nodes wich containts projected vector
IN     global_proj_length    length of (projected) vector
IN     global_proj_width     width of (projected) vector
IN     global_align          alignment to template
IN     templ                 template for vector and matrix distribution
OUT    new_proj_mvector      object describing created projected multi-vector

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_free_enter   (PLA_Obj   *obj)

/*----------------------------------------------------------------------------

Purpose : Linear algebra object destructor.


IN/OUT     object      object to be freed
-----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( *obj != NULL ) 
      if ( (*obj)->pla_obj != NULL ) PLA_Obj_free( &(*obj)->pla_obj );
  }
#endif

  return value;
}


int PLA_Obj_free_exit   (PLA_Obj   *obj)

/*----------------------------------------------------------------------------

Purpose : Linear algebra object destructor.


IN/OUT     obj      object to be freed
-----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_objtype_enter  ( PLA_Obj     obj,    int *objtype)

/*----------------------------------------------------------------------------

Purpose : Extract object type from object.


IN         obj        object to be queried
OUT        objtype    object type of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_objtype";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( objtype == NULL ){
      PLA_Warning( "objtype has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_objtype_exit  ( PLA_Obj     obj,    int *objtype)

/*----------------------------------------------------------------------------

Purpose : Extract object type from object.


IN         obj        object to be queried
OUT        objtype    object type of object

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
/*    int objtype_R12;

    PLA_Obj_objtype( obj->pla_obj, &objtype_R12 );
    if ( *objtype != objtype_R12 ) {
      PLA_Warning( "object type inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Obj_datatype_enter  (PLA_Obj     obj,    MPI_Datatype *datatype)

/*----------------------------------------------------------------------------

Purpose : Extract datatype from object.


IN         obj         object to be queried
OUT        datatype    datatype of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_datatype";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( datatype == NULL ){
      PLA_Warning( "datatype has illegal address" );
      value--;
    }      
  }

  return value;
}

int PLA_Obj_datatype_exit  (PLA_Obj     obj,    MPI_Datatype *datatype)

/*----------------------------------------------------------------------------

Purpose : Extract datatype from object.


IN         obj         object to be queried
OUT        datatype    datatype of object

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
/*    MPI_Datatype datatype_R12;

    PLA_Obj_datatype( obj->pla_obj, &datatype_R12 );
    if ( *datatype != datatype_R12 ) {
      PLA_Warning( "data type inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_template_enter   (PLA_Obj   obj,     PLA_Template    *templ )

/*----------------------------------------------------------------------------

Purpose : Extract template from object.


IN         obj               object to be queried
OUT        templ             template

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_template";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( templ == NULL ){
      PLA_Warning( "templ has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_template_exit   (PLA_Obj   obj,     PLA_Template    *templ )

/*----------------------------------------------------------------------------

Purpose : Extract template from object.


IN         obj               object to be queried
OUT        templ             template

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
/*    PLA_Template template_R12;

    PLA_Obj_template( obj->pla_obj, &template_R12 );
    if ( *templ != template_R12 ) {
      PLA_Warning( "template inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_global_info_enter   (PLA_Obj   obj,             int    *global_length,
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
  char 
    routine_name[ 35 ] = "PLA_Obj_global_info";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
  }

  return value;
}


int PLA_Obj_global_info_exit   (PLA_Obj   obj,             int    *global_length,
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
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_global_length_enter (PLA_Obj   obj,             int    *global_length )

/*----------------------------------------------------------------------------

Purpose :  Extract global length from linear algebra object.

IN         obj                 object to be queried
OUT        global_length       global row dimension of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_global_length";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( global_length == NULL ){
      PLA_Warning( "global_length has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_global_length_exit (PLA_Obj   obj,             int    *global_length )

/*----------------------------------------------------------------------------

Purpose :  Extract global length from linear algebra object.

IN         obj                 object to be queried
OUT        global_length       global row dimension of object

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
/*    int global_length_R12;

    PLA_Obj_global_length( obj->pla_obj, &global_length_R12 );
    if ( *global_length != global_length_R12 ) {
      PLA_Warning( "global_length inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Obj_global_width_enter (PLA_Obj   obj,             int    *global_width )

/*----------------------------------------------------------------------------

Purpose :  Extract global width from linear algebra object.

IN         obj                 object to be queried
OUT        global_width        global column dimension of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_global_width";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( global_width == NULL ){
      PLA_Warning( "global_width has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_global_width_exit (PLA_Obj   obj,             int    *global_width )

/*----------------------------------------------------------------------------

Purpose :  Extract global width from linear algebra object.

IN         obj                 object to be queried
OUT        global_width        global column dimension of object

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
/*    int global_width_R12;

    PLA_Obj_global_width( obj->pla_obj, &global_width_R12 );
    if ( *global_width != global_width_R12 ) {
      PLA_Warning( "global_width inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_project_onto_enter  (PLA_Obj   obj,         int    *project_onto )

/*----------------------------------------------------------------------------

Purpose :  Extract projection information from linear algebra object.


IN         obj                 object to be queried
OUT        project_onto        direction of projection

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_project_onto";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( project_onto == NULL ){
      PLA_Warning( "project_onto has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_project_onto_exit  (PLA_Obj   obj,         int    *project_onto )

/*----------------------------------------------------------------------------

Purpose :  Extract projection information from linear algebra object.


IN         obj                 object to be queried
OUT        project_onto        direction of projection

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
/*    int project_onto_R12;

    PLA_Obj_project_onto( obj->pla_obj, &project_onto_R12 );
    if ( *project_onto != project_onto_R12 ) {
      PLA_Warning( "project_onto inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_owner_row_enter     (PLA_Obj   obj,           int       *owner_row )

/*----------------------------------------------------------------------------

Purpose :  Extract mesh row owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_row           row index of owning node(s)

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_owner_row";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( owner_row == NULL ){
      PLA_Warning( "owner_row has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_owner_row_exit     (PLA_Obj   obj,           int       *owner_row )

/*----------------------------------------------------------------------------

Purpose :  Extract mesh row owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_row           row index of owning node(s)

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_owner_row_consistent_with_R12 ( obj ) ) value--;
/*    int owner_row_R12;

    PLA_Obj_owner_row( obj->pla_obj, &owner_row_R12 );
    if ( *owner_row != owner_row_R12 ) {
      PLA_Warning( "owner_row inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_owner_col_enter     (PLA_Obj   obj,           int       *owner_col )

/*----------------------------------------------------------------------------

Purpose :  Extract mesh column owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_col           column index of owning node(s)

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_owner_col";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( owner_col == NULL ){
      PLA_Warning( "owner_col has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_owner_col_exit     (PLA_Obj   obj,           int       *owner_col )

/*----------------------------------------------------------------------------

Purpose :  Extract mesh column owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_col           column index of owning node(s)

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_owner_col_consistent_with_R12 ( obj ) ) value--;
/*    int owner_col_R12;

    PLA_Obj_owner_col( obj->pla_obj, &owner_col_R12 );
    if ( *owner_col != owner_col_R12 ) {
      PLA_Warning( "owner_col inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_global_align_enter  (PLA_Obj   obj,      int       *global_align )

/*----------------------------------------------------------------------------

WARNING: check this out!!

Purpose :  Extract global alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align        Alignment to template

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_global_align";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( global_align == NULL ){
      PLA_Warning( "global_align has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_global_align_exit  (PLA_Obj   obj,      int       *global_align )

/*----------------------------------------------------------------------------

WARNING: check this out!!

Purpose :  Extract global alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align        Alignment to template

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_global_align_row_consistent_with_R12 ( obj ) ) value--;
/*    int global_align_R12;

    PLA_Obj_global_align( obj->pla_obj, &global_align_R12 );
    if ( *global_align != global_align_R12 ) {
      PLA_Warning( "global_align inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_global_align_row_enter 
                    (PLA_Obj   obj,      int       *global_align_row )

/*-------------------------------------------------------------------------

Purpose :  Extract global row alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_row  row alignment to template

-------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_global_align_row";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( global_align_row == NULL ){
      PLA_Warning( "global_align_row has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_global_align_row_exit       
                    (PLA_Obj   obj,      int       *global_align_row )

/*-------------------------------------------------------------------------

Purpose :  Extract global row alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_row  row alignment to template

-------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_global_align_row_consistent_with_R12 ( obj ) ) value--;
/*    int global_align_row_R12;

    PLA_Obj_global_align_row( obj->pla_obj, &global_align_row_R12 );
    if ( *global_align_row != global_align_row_R12 ) {
      PLA_Warning( "global_align_row inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_global_align_col_enter
          (PLA_Obj   obj,      int       *global_align_col )

/*-------------------------------------------------------------------------

Purpose :  Extract global column alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_col    column alignment to template

-------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_global_align_col";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( global_align_col == NULL ){
      PLA_Warning( "global_align_col has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_global_align_col_exit
          (PLA_Obj   obj,      int       *global_align_col )

/*-------------------------------------------------------------------------

Purpose :  Extract global column alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_col    column alignment to template

-------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_global_align_col_consistent_with_R12 ( obj ) ) value--;
/*    int global_align_col_R12;

    PLA_Obj_global_align_col( obj->pla_obj, &global_align_col_R12 );
    if ( *global_align_col != global_align_col_R12 ) {
      PLA_Warning( "global_align_col insistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_local_info_enter   (PLA_Obj   obj,             int    *local_length,
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
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Obj_local_info_exit   (PLA_Obj   obj,             int    *local_length,
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
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_local_length_enter   (
	    PLA_Obj   obj,       int    *local_length )

/*----------------------------------------------------------------------------

Purpose :  Extract local length from linear algebra object.

IN         obj                 object to be queried
OUT        local_length        local row dimension of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_local_length";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( local_length == NULL ){
      PLA_Warning( "local_length has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_local_length_exit   (
	    PLA_Obj   obj,       int    *local_length )

/*----------------------------------------------------------------------------

Purpose :  Extract local length from linear algebra object.

IN         obj                 object to be queried
OUT        local_length        local row dimension of object

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_local_length_consistent_with_R12 ( obj ) ) value--;
/*
    int local_length_R12;

    PLA_Obj_local_length( obj->pla_obj, &local_length_R12 );
    if ( *local_length != local_length_R12 ) {
      if ( *local_length != 0 || local_length_R12 != PLA_NO_DIMENSION ) {
	PLA_Warning( "local_length inconsistent between R2.0 and R1.2" );
	printf("%d vs %d\n", *local_length, local_length_R12 );
	value--;
      }
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_local_width_enter   (PLA_Obj   obj,    int       *local_width )

/*----------------------------------------------------------------------------

Purpose :  Extract local width from linear algebra object.

IN         obj                 object to be queried
OUT        local_width         local column dimension of object

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_local_width";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( local_width == NULL ){
      PLA_Warning( "local_width has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_local_width_exit   (PLA_Obj   obj,    int       *local_width )

/*----------------------------------------------------------------------------

Purpose :  Extract local width from linear algebra object.

IN         obj                 object to be queried
OUT        local_width         local column dimension of object

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_local_width_consistent_with_R12 ( obj ) ) value--;
/*    int local_width_R12;

    PLA_Obj_local_width( obj->pla_obj, &local_width_R12 );
    if ( *local_width != local_width_R12 ) {
      if ( *local_width != 0 || local_width_R12 != PLA_NO_DIMENSION ) {
	PLA_Warning( "local_width inconsistent between R2.0 and R1.2" );
	printf("%d vs %d\n", *local_width, local_width_R12 );
	value--;
      }
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_local_buffer_enter  (PLA_Obj   obj,    void   **local_buffer )

/*----------------------------------------------------------------------------

Purpose :  Extract local buffer from linear algebra object.

IN         obj                 object to be queried
OUT        local_buffer        address of local data

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_local_buffer";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( local_buffer == NULL ){
      PLA_Warning( "local_buffer has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_local_buffer_exit  (PLA_Obj   obj,    void   **local_buffer )

/*----------------------------------------------------------------------------

Purpose :  Extract local buffer from linear algebra object.

IN         obj                 object to be queried
OUT        local_buffer        address of local data

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_local_ldim_enter   (PLA_Obj   obj,             int    *local_ldim )

/*----------------------------------------------------------------------------

Purpose :  Extract local leading dimension from linear algebra object.

IN         obj                 object to be queried
OUT        local_ldim          leading dimension of array holding local data

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_local_ldim";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( local_ldim == NULL ){
      PLA_Warning( "local_ldim has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_local_ldim_exit   (PLA_Obj   obj,             int    *local_ldim )

/*----------------------------------------------------------------------------

Purpose :  Extract local leading dimension from linear algebra object.

IN         obj                 object to be queried
OUT        local_ldim          leading dimension of array holding local data

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_local_ldim_consistent_with_R12 ( obj ) ) value--;
/*    int local_ldim_R12;

    PLA_Obj_local_ldim( obj->pla_obj, &local_ldim_R12 );
    if ( *local_ldim != local_ldim_R12 ) {
      PLA_Warning( "local_ldim inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_local_stride_enter   (PLA_Obj   obj,        int    *local_stride )

/*----------------------------------------------------------------------------

Purpose :  Extract local stride from linear algebra object.

IN         obj                 object to be queried
OUT        local_stride        stride between entries in a column or vector

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_local_stride";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( local_stride == NULL ){
      PLA_Warning( "local_stride has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_local_stride_exit   (PLA_Obj   obj,        int    *local_stride )

/*----------------------------------------------------------------------------

Purpose :  Extract local stride from linear algebra object.

IN         obj                 object to be queried
OUT        local_stride        stride between entries in a column or vector

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_local_stride_consistent_with_R12 ( obj ) ) value--;
/*    int local_stride_R12;

    PLA_Obj_local_stride( obj->pla_obj, &local_stride_R12 );
    if ( *local_stride != local_stride_R12 ) {
      PLA_Warning( "local_stride inconsistent between R2.0 and R1.2" );
      value--;
    } */
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_get_local_contents_enter   (
        PLA_Obj      obj,              int      trans,
        int           *rows_in_buf,     int      *cols_in_buf,
        void          *buf,             int      leading_dim_buf,
        int           stride_buf)

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
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_get_local_contents";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( trans != PLA_NO_TRANS && trans != PLA_TRANS ){
      PLA_Warning( "trans has illegal value" );
      value--;
    }      

    if ( rows_in_buf == NULL ) {
      PLA_Abort( "rows_in_buf has illegal address (NULL)", __LINE__, __FILE__ );
      value--;
    }      

    if ( cols_in_buf == NULL ) {
      PLA_Abort( "cols_in_buf has illegal address (NULL)", __LINE__, __FILE__ );
      value--;
    }      


    if ( buf == NULL ) {
      int length, width;

      PLA_OBJ_LOCAL_LENGTH( obj, &length );
      PLA_OBJ_LOCAL_WIDTH( obj, &width );
      if ( length != 0 && width != 0 ) {
	PLA_Abort( "buf has illegal address (NULL)", __LINE__, __FILE__ );
	value--;
      }
    }      

    if ( trans == PLA_NO_TRANS ){
      int length;

      PLA_OBJ_LOCAL_LENGTH( obj, &length );
      if ( length > leading_dim_buf ) {
	PLA_Warning( "leading_dim_buf less than row dimension" );
	value--;
      }
    }      
    else if ( trans == PLA_TRANS ){
      int width;

      PLA_OBJ_LOCAL_WIDTH( obj, &width );
      if ( width > leading_dim_buf ) {
	PLA_Warning( "leading_dim_buf less than column dimension" );
	value--;
      }
    }      
    
    if ( stride_buf != 1 ){
      PLA_Warning( "stride_buf != 1 not implemented" );
      value--;
    }
  }

  return value;
}


int PLA_Obj_get_local_contents_exit   (
             PLA_Obj       obj,             int      trans,
             int           *rows_in_buf,     int      *cols_in_buf,
             void          *buf,             int      leading_dim_buf,
             int           stride_buf )

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
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_set_local_contents_enter   (
              int      trans,                 int      rows_in_buf,
              int      cols_in_buf,           void     *buf,
              int      leading_dim_buf,       int      stride_buf,
              PLA_Obj  obj )

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
    value = PLA_SUCCESS, length, width;
  char 
    routine_name[ 35 ] = "PLA_Obj_set_local_contents";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( trans != PLA_NO_TRANS && trans != PLA_TRANS ){
      PLA_Warning( "trans has illegal value" );
      value--;
    }      

    if ( rows_in_buf < 0 ) {
      PLA_Warning( "negative rows_in_buf" );
      value--;
    }      

    if ( cols_in_buf < 0 ) {
      PLA_Warning( "negative cols_in_buf" );
      value--;
    }      

    PLA_OBJ_LOCAL_LENGTH( obj, &length );
    PLA_OBJ_LOCAL_WIDTH( obj, &width );

    if ( buf == NULL ) {
      if ( length != 0 && width != 0 ) {
	PLA_Abort( "buf has illegal address (NULL)", __LINE__, __FILE__ );
	value--;
      }
    }      


    if ( trans == PLA_NO_TRANS ){
      if ( length != rows_in_buf ) {
	PLA_Warning( "row dimensions of obj and buffer don't match" );
	value--;
      }

      if ( width != cols_in_buf ) {
	PLA_Warning( "column dimensions of obj and buffer don't match" );
	value--;
      }

      if ( length > leading_dim_buf ) {
	PLA_Warning( "leading_dim_buf less than row dimension" );
	value--;
      }
    }      
    else if ( trans == PLA_TRANS ){
      if ( width != rows_in_buf ) {
	PLA_Warning( "dimensions of obj and buffer don't match" );
	value--;
      }

      if ( length != cols_in_buf ) {
	PLA_Warning( "dimensions of obj and buffer don't match" );
	value--;
      }

      if ( width > leading_dim_buf ) {
	PLA_Warning( "leading_dim_buf less than column dimension" );
	value--;
      }
    }      
    
    if ( stride_buf != 1 ){
      PLA_Warning( "stride_buf != 1 not implemented" );
      value--;
    }
  }

  return value;
}


int PLA_Obj_set_local_contents_exit   (
           int      trans,                 int      rows_in_buf,
           int      cols_in_buf,           void     *buf,
           int      leading_dim_buf,       int      stride_buf,
           PLA_Obj  obj )

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
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Obj_set_local_contents( 
         ( trans == PLA_TRANS ? PLA_TRANS : PLA_NO_TRANS ), 
         rows_in_buf, cols_in_buf, buf, 
         leading_dim_buf, stride_buf, obj->pla_obj );

    if ( !PLA_Obj_contents_consistent_with_R12( obj ) ){
      PLA_Warning( "contents inconsistent between R2.0 and R1.2" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Obj_set_enter   (PLA_Obj      obj,    MPI_Datatype   datatype,
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
    return_value = PLA_SUCCESS, length, width;
  char 
    routine_name[ 35 ] = "PLA_Obj_set";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      return_value--;
    }
    
    if ( value == NULL ) {
      PLA_Abort( "value has illegal address (NULL)", __LINE__, __FILE__ );
      return_value--;
    }      
  }

  return return_value;
}


int PLA_Obj_set_exit   (PLA_Obj      obj,    MPI_Datatype   datatype,
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
  int return_value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Obj_set( obj->pla_obj, datatype, value );

    if ( !PLA_Obj_contents_consistent_with_R12( obj ) ){
      PLA_Warning( "contents inconsistent between R2.0 and R1.2" );
      return_value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return return_value;
}


/***************************************************************************/

int PLA_Obj_set_to_zero_enter   (PLA_Obj      obj)

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to zero (0).
The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE, etc.)
where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/
{
  int
    return_value = PLA_SUCCESS, length, width;
  char 
    routine_name[ 35 ] = "PLA_Obj_set_to_zero";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      return_value--;
    }
  }

  return return_value;
}


int PLA_Obj_set_to_zero_exit   (PLA_Obj      obj)

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to zero (0).
The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE, etc.)
where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/
{
  int return_value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Obj_set_to_zero( obj->pla_obj );

    if ( !PLA_Obj_contents_consistent_with_R12( obj ) ){
      PLA_Warning( "contents inconsistent between R2.0 and R1.2" );
      return_value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return return_value;
}


/***************************************************************************/

int PLA_Obj_set_to_one_enter   (PLA_Obj      obj)

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to one (1).
The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE, etc.)
where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/
{
  int
    return_value = PLA_SUCCESS, length, width;
  char 
    routine_name[ 35 ] = "PLA_Obj_set_to_one";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      return_value--;
    }
  }

  return return_value;
}


int PLA_Obj_set_to_one_exit   (PLA_Obj      obj)

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to one (1).
The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE, etc.)
where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/
{
  int return_value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Obj_set_to_one( obj->pla_obj );

    if ( !PLA_Obj_contents_consistent_with_R12( obj ) ){
      PLA_Warning( "contents inconsistent between R2.0 and R1.2" );
      return_value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return return_value;
}


/***************************************************************************/

int PLA_Obj_set_to_minus_one_enter   (PLA_Obj      obj)

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to negative
one (-1). The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE,
etc.) where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/
{
  int
    return_value = PLA_SUCCESS, length, width;
  char 
    routine_name[ 35 ] = "PLA_Obj_set_to_minus_one";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      return_value--;
    }
  }

  return return_value;
}


int PLA_Obj_set_to_minus_one_exit   (PLA_Obj      obj)

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to negative
one (-1). The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE,
etc.) where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/
{
  int return_value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Obj_set_to_minus_one( obj->pla_obj );

    if ( !PLA_Obj_contents_consistent_with_R12( obj ) ){
      PLA_Warning( "contents inconsistent between R2.0 and R1.2" );
      return_value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

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

int PLA_Obj_view_enter   (PLA_Obj   old_obj,       int      global_length,
                    int         global_width,    int      align_row,
                    int         align_col,       PLA_Obj *new_obj)

/*----------------------------------------------------------------------------

Purpose : Create a view (reference) into an existing linear algebra
object.


IN     old_obj          object into which view is taken
IN     global_length    row dimension of view
IN     global_width     column dimension of view
IN     align_row        row index in old object of upper-left-hand element of view
IN     align_col        column index in old object of upper-left-hand element of view
IN/OUT new_obj          created view

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Obj_view_exit   (PLA_Obj   old_obj,       int      global_length,
                    int         global_width,    int      align_row,
                    int         align_col,       PLA_Obj *new_obj)

/*----------------------------------------------------------------------------

Purpose : Create a view (reference) into an existing linear algebra
object.


IN     old_obj          object into which view is taken
IN     global_length    row dimension of view
IN     global_width     column dimension of view
IN     align_row        row index in old object of upper-left-hand element of view
IN     align_col        column index in old object of upper-left-hand element of view
IN/OUT new_obj          created view

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_view_all_enter   (PLA_Obj     old_obj,      PLA_Obj     *new_obj)

/*----------------------------------------------------------------------------

Purpose : Create a view (reference) into all of existing linear algebra
object.


IN          old_obj          object into which view is taken
IN/OUT      new_obj          created view

----------------------------------------------------------------------------*/
{
  int
    return_value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_view_all";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( old_obj == NULL || !PLA_Valid_object( old_obj ) ) {
      PLA_Warning( "Invalid object" );
      return_value--;
    }

    if ( *new_obj != NULL && !PLA_Valid_object( *new_obj ) ) {
      PLA_Warning( "Invalid object" );
      return_value--;
    }
  }

  return return_value;
}



int PLA_Obj_view_all_exit   (PLA_Obj     old_obj,      PLA_Obj     *new_obj)

/*----------------------------------------------------------------------------

Purpose : Create a view (reference) into all of existing linear algebra
object.


IN          old_obj          object into which view is taken
IN/OUT      new_obj          created view

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( *new_obj != old_obj ){
      (*new_obj)->pla_obj = NULL;
      PLA_Obj_view_all( old_obj->pla_obj, &(*new_obj)->pla_obj );
    }

    if ( !PLA_Obj_consistent_with_R12( *new_obj ) ){
      PLA_Warning( "Objects inconsistent" );
      value--;
    }

    if ( !PLA_Obj_contents_consistent_with_R12( *new_obj ) ){
      PLA_Warning( "Object contents inconsistent" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_view_swap_enter    ( PLA_Obj     *obj1,        PLA_Obj   *obj2)

/*----------------------------------------------------------------------------

Purpose : Swaps two views.


IN/OUT      obj1          linear algebra object
IN/OUT      obj2          linear algebra object

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Obj_view_swap_exit    ( PLA_Obj     *obj1,        PLA_Obj   *obj2)

/*----------------------------------------------------------------------------

Purpose : Swaps two views.


IN/OUT      obj1          linear algebra object
IN/OUT      obj2          linear algebra object

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_split_4_enter  (
     PLA_Obj       obj,
     int           length,            int           width,
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
    count_equal_obj = 0,
    global_length, global_width;
  
  char 
    routine_name[ 35 ] = "PLA_Obj_split_4";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      return_value--;
    }

    PLA_OBJ_GLOBAL_LENGTH( obj, &global_length );
    if ( ( length >= 0 ? length : -length ) > global_length ){
      PLA_Warning( "row dimension of split larger than global length" );
      return_value--;
    }

    PLA_OBJ_GLOBAL_WIDTH( obj, &global_width );
    if ( ( width >= 0 ? width : -width ) > global_width ){
      PLA_Warning( "column dimension of split larger than global width" );
      return_value--;
    }
      
    if ( upper_left_obj != PLA_DUMMY &&
         *upper_left_obj != NULL && !PLA_Valid_object( *upper_left_obj ) ) {
      PLA_Warning( "Invalid upper_left object" );
      return_value--;
    }

    if ( upper_right_obj != PLA_DUMMY &&
         *upper_right_obj != NULL && !PLA_Valid_object( *upper_right_obj ) ) {
      PLA_Warning( "Invalid upper_right object" );
      return_value--;
    }

    if ( lower_left_obj != PLA_DUMMY &&
         *lower_left_obj != NULL && !PLA_Valid_object( *lower_left_obj ) ) {
      PLA_Warning( "Invalid lower_left object" );
      return_value--;
    }

    if ( lower_right_obj != PLA_DUMMY &&
         *lower_right_obj != NULL && !PLA_Valid_object( *lower_right_obj ) ) {
      PLA_Warning( "Invalid lower_right object" );
      return_value--;
    }

    if ( upper_right_obj != PLA_DUMMY && *upper_right_obj == obj ) 
      count_equal_obj++;
    if ( upper_left_obj != PLA_DUMMY && *upper_left_obj == obj )  count_equal_obj++;
    if ( lower_right_obj != PLA_DUMMY && *lower_right_obj == obj ) count_equal_obj++;
    if ( lower_left_obj != PLA_DUMMY && *lower_left_obj == obj )  count_equal_obj++;
    if ( count_equal_obj > 1 ){
      PLA_Warning( "more than one output object equal to input object" );
      return_value--;
    }
  }

  return return_value;
}


int PLA_Obj_split_4_exit  (
     PLA_Obj       obj,
     int           length,            int           width,
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
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( upper_left_obj != PLA_DUMMY ) {
      PLA_Obj_split_4( (*upper_left_obj)->pla_obj, length, width,
		      &(*upper_left_obj)->pla_obj, PLA_DUMMY,
                      PLA_DUMMY,                   PLA_DUMMY );

      if ( !PLA_Obj_consistent_with_R12( *upper_left_obj ) ){
	PLA_Warning( "Upper_left object inconsistent" );
	value--;
      }

      if ( !PLA_Obj_contents_consistent_with_R12( *upper_left_obj ) ){
	PLA_Warning( "contents of upper_left_obj inconsistent between R2.0 and R1.2" );
	value--;
      }
    }

    if ( upper_right_obj != PLA_DUMMY ) {
      PLA_Obj_split_4( (*upper_right_obj)->pla_obj, length, width,
		      PLA_DUMMY,     &(*upper_right_obj)->pla_obj,
                      PLA_DUMMY,     PLA_DUMMY );

      if ( !PLA_Obj_consistent_with_R12( *upper_right_obj ) ){
	PLA_Warning( "Upper_right object inconsistent" );
	value--;
      }

      if ( !PLA_Obj_contents_consistent_with_R12( *upper_right_obj ) ){
	PLA_Warning( "contents of upper_right_obj inconsistent between R2.0 and R1.2" );
	value--;
      }
    }

    if ( lower_left_obj != PLA_DUMMY ) {
      PLA_Obj_split_4( (*lower_left_obj)->pla_obj, length, width,
                      PLA_DUMMY,                   PLA_DUMMY,
		      &(*lower_left_obj)->pla_obj, PLA_DUMMY);

      if ( !PLA_Obj_consistent_with_R12( *lower_left_obj ) ){
	PLA_Warning( "lower_left object inconsistent" );
	value--;
      }

      if ( !PLA_Obj_contents_consistent_with_R12( *lower_left_obj ) ){
	PLA_Warning( "contents of lower_left_obj inconsistent between R2.0 and R1.2" );
	value--;
      }
    }

    if ( lower_right_obj != PLA_DUMMY ) {
      PLA_Obj_split_4( (*lower_right_obj)->pla_obj, length, width,
                      PLA_DUMMY,     PLA_DUMMY,
		      PLA_DUMMY,     &(*lower_right_obj)->pla_obj );

      if ( !PLA_Obj_consistent_with_R12( *lower_right_obj ) ){
	PLA_Warning( "lower_right object inconsistent" );
	value--;
      }
      if ( !PLA_Obj_contents_consistent_with_R12( *lower_right_obj ) ){
	PLA_Warning( "contents of lower_right_obj inconsistent between R2.0 and R1.2" );
	value--;
      }
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}


/***************************************************************************/

int PLA_Obj_horz_split_2_enter   (  PLA_Obj     obj,      int          length,
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
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Obj_horz_split_2_exit   (  PLA_Obj     obj,      int          length,
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
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_vert_split_2_enter   (  PLA_Obj     obj,           int          width,
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
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Obj_vert_split_2_exit   (  PLA_Obj     obj,           int          width,
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
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_view_shift_enter   (PLA_Obj      obj,           int      length_top,
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
  int return_value = PLA_SUCCESS;

  char 
    routine_name[ 35 ] = "PLA_Obj_view_shift";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      return_value--;
    }
  }

  return return_value;
}


int PLA_Obj_view_shift_exit   (PLA_Obj      obj,           int      length_top,
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
  int value = PLA_SUCCESS;

  char 
    routine_name[ 35];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Obj_view_shift( obj->pla_obj, length_top, width_left, width_right, 
                        length_bottom );

    if ( !PLA_Obj_consistent_with_R12( obj ) ){
      PLA_Warning( "View shift of object inconsistent" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Obj_split_size_enter   (PLA_Obj      obj,      int      side,
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
  char 
    routine_name[ 35 ] = "PLA_Obj_split_size";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( size == NULL ){
      PLA_Warning( "size has illegal address" );
      value--;
    }      

    if ( owner == NULL ){
      PLA_Warning( "owner has illegal address" );
      value--;
    }      
  }

  return value;
}


int PLA_Obj_split_size_exit   (PLA_Obj      obj,      int      side,
                          int          *size,    int      *owner)

/*----------------------------------------------------------------------------

Purpose : Compute size of split to next block boundary.

IN       obj       object to be split
IN       side      side of split
OUT      size      size to next template subblock split
OUT      owner     index of row or column of nodes that owns the split block

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( !PLA_Obj_sizes_and_owners_consistent_with_R12 ( obj ) ) value--;
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Vector_create_conf_to_enter   (PLA_Obj      obj,        PLA_Obj     *new_vector)

/*----------------------------------------------------------------------------

Purpose : Create vector conformal to given object.


IN      obj          original object
OUT     new_vector   created object

-------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Vector_create_conf_to_exit   (PLA_Obj      obj,        PLA_Obj     *new_vector)

/*----------------------------------------------------------------------------

Purpose : Create vector conformal to given object.


IN      obj          original object
OUT     new_vector   created object

-------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Mvector_create_conf_to_enter   (PLA_Obj      obj,        int      global_width,
                                  PLA_Obj     *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create mvector conformal to given object.


IN      obj                original object
IN      global_width       number of vectors in mvector
OUT     new_obj            created object

-------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS,
    objtype;
  char 
    routine_name[ 35 ] = "PLA_Mvector_create_conf_to";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid input object parameter" );
      value--;
    }

    PLA_OBJ_OBJTYPE( obj, &objtype );
    if ( objtype != PLA_MVECTOR && objtype != PLA_MATRIX &&
	 objtype != PLA_PMVECTOR ){
      PLA_Warning( "Cannot create Mvector conf. to given objtype" );
      value--;
    }

    if ( global_width < 0 ){
      PLA_Warning( "Negative global width" );
      value--;
    }

    if ( *new_obj != NULL && !PLA_Valid_object( *new_obj ) ) {
      PLA_Warning( "Invalid output object parameter" );
      value--;
    }
  }

  return value;
}


int PLA_Mvector_create_conf_to_exit   (PLA_Obj      obj,        int      global_width,
                                  PLA_Obj     *new_obj )

/*----------------------------------------------------------------------------

Purpose : Create mvector conformal to given object.


IN      obj                original object
IN      global_width       number of vectors in mvector
OUT     new_obj            created object

-------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;

  char 
    routine_name[ 35];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    PLA_Mvector_create_conf_to( 
		     obj->pla_obj, global_width, &(*new_obj)->pla_obj );

    if ( !PLA_Obj_consistent_with_R12( *new_obj ) ){
      PLA_Warning( "Objects inconsistent" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Pvector_create_conf_to_enter   (PLA_Obj      obj,        int         project_onto,
                                  int          owner,      PLA_Obj     *new_pvector)

/*----------------------------------------------------------------------------

Purpose : Create projected vector conformal to given object.


IN      obj                original object
IN      project_onto       mesh direction onto which to project
IN      owner              index of row or column of nodes which contains owning nodes
OUT     new_pvector        created object

-------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Pvector_create_conf_to_exit   (PLA_Obj      obj,        int         project_onto,
                                  int          owner,      PLA_Obj     *new_pvector)

/*----------------------------------------------------------------------------

Purpose : Create projected vector conformal to given object.


IN      obj                original object
IN      project_onto       mesh direction onto which to project
IN      owner              index of row or column of nodes which contains owning nodes
OUT     new_pvector        created object

-------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Pmvector_create_conf_to_enter   (PLA_Obj      obj,        int         project_onto,
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
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Pmvector_create_conf_to_exit   (PLA_Obj      obj,        int         project_onto,
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
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Matrix_create_conf_to_enter   (PLA_Obj      matrix,        PLA_Obj     *new_matrix)

/*----------------------------------------------------------------------------

Purpose : Create matrix conformal to given object (must be a matrix).



IN      matrix       original object
OUT     new_matrix   created object

-------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Matrix_create_conf_to_exit   (PLA_Obj      matrix,        PLA_Obj     *new_matrix)

/*----------------------------------------------------------------------------

Purpose : Create matrix conformal to given object (must be a matrix).



IN      matrix       original object
OUT     new_matrix   created object

-------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Mscalar_create_conf_to_enter   (PLA_Obj      obj,        int         owner_row,
                                  int          owner_col,  PLA_Obj     *new_mscalar)

/*----------------------------------------------------------------------------

Purpose : Create multiscalar conformal to given object.


IN      obj           original object
IN      owner_row     index of row of nodes that contains owning node(s)
IN      owner_col     index of column of nodes that contains owning node(s)
OUT     new_mscalar   created object (multiscalar)

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}


int PLA_Mscalar_create_conf_to_exit   (PLA_Obj      obj,        int         owner_row,
                                  int          owner_col,  PLA_Obj     *new_mscalar)

/*----------------------------------------------------------------------------

Purpose : Create multiscalar conformal to given object.


IN      obj           original object
IN      owner_row     index of row of nodes that contains owning node(s)
IN      owner_col     index of column of nodes that contains owning node(s)
OUT     new_mscalar   created object (multiscalar)

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;

  return value;
}

/***************************************************************************/

int PLA_Obj_objtype_cast_enter   (PLA_Obj      obj,         int      obj_type)

/*----------------------------------------------------------------------------

Purpose : Cast the linear algebra object to have given object type.


IN/OUT        obj            original object
IN            objtype        new object type

-------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_objtype_cast";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid input object parameter" );
      value--;
    }

    if ( obj_type != PLA_MVECTOR && obj_type != PLA_MATRIX &&
	 obj_type != PLA_PMVECTOR && obj_type != PLA_MSCALAR ){
      PLA_Warning( "Illegal objtype" );
      value--;
    }
  }

  return value;
}


int PLA_Obj_objtype_cast_exit   (PLA_Obj      obj,         int      obj_type)

/*----------------------------------------------------------------------------

Purpose : Cast the linear algebra object to have given object type.


IN/OUT        obj            original object
IN            obj_type        new object type

-------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;

  char 
    routine_name[ 35];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    switch ( obj_type ) {
    case PLA_MSCALAR:
      PLA_Obj_objtype_cast( obj->pla_obj, PLA_MSCALAR );
      break;
    case PLA_MATRIX:
      PLA_Obj_objtype_cast( obj->pla_obj, PLA_MATRIX );
      break;
    case PLA_MVECTOR:
      PLA_Obj_objtype_cast( obj->pla_obj, PLA_MVECTOR );
      break;
    case PLA_PMVECTOR:
      PLA_Obj_objtype_cast( obj->pla_obj, PLA_PMVECTOR );
      break;
    }

    if ( !PLA_Obj_consistent_with_R12( obj ) ){
      PLA_Warning( "Objects inconsistent" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Obj_set_orientation_enter   (PLA_Obj     obj,     int     project_onto)

/*----------------------------------------------------------------------------

Purpose : Set the orientation of a matrix panel.


IN/OUT      obj              original object
IN          project_onto     annotation indicating orientation

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_set_orientation";

  PLA_Routine_stack_push( routine_name );

#if DEBUG==1
  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( project_onto != PLA_PROJ_ONTO_ROW && 
         project_onto != PLA_PROJ_ONTO_COL ) {
      PLA_Warning( "project_onto has illegal value" );
      value--;
    }
  }
#endif

  return value;
}


int PLA_Obj_set_orientation_exit   (PLA_Obj     obj,     int     project_onto)

/*----------------------------------------------------------------------------

Purpose : Set the orientation of a matrix panel.


IN/OUT      obj              original object
IN          project_onto     annotation indicating orientation

----------------------------------------------------------------------------*/
{
  int 
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

#if DEBUG==1
  if ( PLA_CHECK_AGAINST_R12 ){
    if ( project_onto == PLA_PROJ_ONTO_COL )
      PLA_Obj_set_orientation( obj->pla_obj, PLA_PROJ_ONTO_COL );
    else
      PLA_Obj_set_orientation( obj->pla_obj, PLA_PROJ_ONTO_ROW );

    if ( !PLA_Obj_consistent_with_R12( obj ) ){
      PLA_Warning( "Objects inconsistent" );
      value--;
    }
  }
#endif

  PLA_Routine_stack_pop( routine_name );

  return value;
}

/***************************************************************************/

int PLA_Obj_get_orientation_enter   (PLA_Obj     obj,     int     *project_onto)


/*----------------------------------------------------------------------------

Purpose : Return the orientation of a matrix panel.


IN/OUT      obj              original object
IN          project_onto     annotation indicating orientation

----------------------------------------------------------------------------*/
{
  int
    value = PLA_SUCCESS;
  char 
    routine_name[ 35 ] = "PLA_Obj_get_orientation";

  PLA_Routine_stack_push( routine_name );

  if ( PLA_CHECK_PARAMETERS ){
    if ( obj == NULL || !PLA_Valid_object( obj ) ) {
      PLA_Warning( "Invalid object" );
      value--;
    }
    
    if ( project_onto == NULL ){
      PLA_Abort( "project_onto has illegal address (NULL)", __LINE__, __FILE__ );
      value--;
    }
  }

  return value;
}


int PLA_Obj_get_orientation_exit   (PLA_Obj     obj,     int     *project_onto)


/*----------------------------------------------------------------------------

Purpose : Return the orientation of a matrix panel.


IN/OUT      obj              original object
IN          project_onto     annotation indicating orientation

----------------------------------------------------------------------------*/
{
  int value = PLA_SUCCESS;
  char 
    routine_name[ 35 ];

  PLA_Routine_stack_pop( routine_name );


  return value;
}

