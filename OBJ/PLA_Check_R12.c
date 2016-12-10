/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

#if DEBUG==1
/***************************************************************************/

int PLA_Obj_consistent_with_R12( PLA_Obj obj )
     
/*----------------------------------------------------------------------------

Purpose : Check if object values are consistent with PLAPACK R12

IN     obj        object to be checked

----------------------------------------------------------------------------*/
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;
  MPI_Datatype 
    datatype_R20, datatype_R12;
  PLA_Template
    templ_R20, templ_R12;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  if ( !PLA_Obj_objtype_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_datatype_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_global_length_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_global_width_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_global_align_row_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_global_align_col_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_owner_row_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_owner_col_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_project_onto_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_local_length_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_local_width_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_local_ldim_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_local_stride_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_template_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !PLA_Obj_sizes_and_owners_consistent_with_R12( obj ) ) 
    return_value = FALSE;

  if ( !return_value )
    PLA_Print_stack_history( );

  return return_value;
}

/***************************************************************************/

int PLA_Obj_contents_consistent_with_R12( PLA_Obj obj )
     
/*----------------------------------------------------------------------------

Purpose : Check if object contents are consistent with PLAPACK R12

IN     obj        object to be checked

----------------------------------------------------------------------------*/
{
  int 
    return_value = TRUE,
    typesize,
    ldim_R12, ldim_R20,
    length, width;
  MPI_Datatype 
    datatype;
  PLA_Obj
    pla_obj;
  char *buffer_R20, *buffer_R12;

  /* for the moment, we don't check! */
  return return_value;

  PLA_Obj_datatype( obj, &datatype );
  MPI_Type_size( datatype, &typesize );

#if DEBUG==1
  pla_obj = obj->pla_obj;
#endif

  PLA_Obj_local_ldim( obj,     &ldim_R20 );
  PLA_Obj_local_ldim ( pla_obj, &ldim_R12 );

  PLA_Obj_local_length( obj,     &length );
  PLA_Obj_local_width ( obj,     &width );

  PLA_Obj_local_buffer( obj,     (void **) &buffer_R20 );
  PLA_Obj_local_buffer ( pla_obj, (void **) &buffer_R12 );

  ldim_R20 *= typesize;
  ldim_R12 *= typesize;
  length   *= typesize;

  return_value = PLA_Compare_local_arrays( 
              length, width, buffer_R20, ldim_R20, buffer_R12, ldim_R12 );

  return return_value;
}


/***************************************************************************/

int PLA_Obj_objtype_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_OBJTYPE            ( obj, &value_R20 );
  PLA_Obj_objtype             ( obj->pla_obj, &value_R12 );

  if ( value_R20 == PLA_MVECTOR && value_R12 == PLA_MVECTOR ) return TRUE;
  if ( value_R20 == PLA_MATRIX && value_R12 == PLA_MATRIX ) return TRUE;
  if ( value_R20 == PLA_MSCALAR && value_R12 == PLA_MSCALAR ) return TRUE;
  if ( value_R20 == PLA_PMVECTOR && value_R12 == PLA_PMVECTOR ) return TRUE;

  PLA_Warning( "inconsistent objtype R20 vs R12" );
  printf("%3d: %d vs %d\n", me, value_R20, value_R12 );

  return FALSE;
}

int PLA_Obj_datatype_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    me;

  MPI_Datatype
    datatype_R20, datatype_R12;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_DATATYPE           ( obj, &datatype_R20 );
  PLA_Obj_datatype            ( obj->pla_obj, &datatype_R12 );
  if ( datatype_R20 != datatype_R12 ) {
    PLA_Warning( "inconsistent datatype R20 vs R12" );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_global_length_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_GLOBAL_LENGTH      ( obj, &value_R20 );
  PLA_Obj_global_length       ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 ) {
    PLA_Warning( "inconsistent global_length R20 vs R12" );
    printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_global_width_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_GLOBAL_WIDTH      ( obj, &value_R20 );
  PLA_Obj_global_width       ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 ) {
    PLA_Warning( "inconsistent global_width R20 vs R12" );
    printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_global_align_row_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_GLOBAL_ALIGN_ROW   ( obj, &value_R20 );
  PLA_Obj_global_align_row    ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 ) {
    PLA_Warning( "inconsistent global_align_row R20 vs R12" );
    printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_global_align_col_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_GLOBAL_ALIGN_COL   ( obj, &value_R20 );
  PLA_Obj_global_align_col    ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 ) {
    PLA_Warning( "inconsistent global_align_col R20 vs R12" );
    printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_owner_row_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_OWNER_ROW          ( obj, &value_R20 );
  PLA_Obj_owner_row           ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 ) {
    if ( value_R20 == PLA_ALL_ROWS && 
         ( value_R12 == PLA_ALL_ROWS || value_R12 == PLA_UNDEFINED ) )
      return TRUE;
    PLA_Warning( "inconsistent owner_row R20 vs R12" );
    printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_owner_col_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_OWNER_COL          ( obj, &value_R20 );
  PLA_Obj_owner_col           ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 ) {
    if ( value_R20 == PLA_ALL_COLS && 
         ( value_R12 == PLA_ALL_COLS || value_R12 == PLA_UNDEFINED ) )
      return TRUE;
    PLA_Warning( "inconsistent owner_col R20 vs R12" );
    printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_project_onto_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;
  static int 
    first_time = TRUE;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_PROJECT_ONTO          ( obj, &value_R20 );
  PLA_Obj_project_onto           ( obj->pla_obj, &value_R12 );

  if ( value_R20 == PLA_PROJ_ONTO_ROW && 
       ( value_R12 != PLA_PROJ_ONTO_COL ) )
    return TRUE;
  if ( value_R20 == PLA_PROJ_ONTO_COL &&
       ( value_R12 != PLA_PROJ_ONTO_ROW ) )
    return TRUE;

  if ( first_time ){
    PLA_Warning( "inconsistent project_onto R20 vs R12" );

    PLA_Print_objtype( obj );
    printf("\n");

    if ( value_R20 == PLA_PROJ_ONTO_COL )
      printf("R20: PLA_PROJ_ONTO_COL\n");
    else if ( value_R20 == PLA_PROJ_ONTO_ROW )
      printf("R20: PLA_PROJ_ONTO_ROW\n");
    else
      printf("R20: PLA_PROJ_ONTO_ unknown\n");

    if ( value_R12 == PLA_PROJ_ONTO_COL )
      printf("R12: PLA_PROJ_ONTO_COL\n");
    else if ( value_R12 == PLA_PROJ_ONTO_ROW )
      printf("R12: PLA_PROJ_ONTO_ROW\n");
    else
      printf("R12: PLA_PROJ_ONTO_ unknown\n");
    
    printf("since this may be due to an inherent inconsistency\n");
    printf("this message is printed only once\n");
           
    first_time = FALSE;

    return FALSE;
  }
  else
    return TRUE;
}

int PLA_Obj_local_length_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_LOCAL_LENGTH      ( obj, &value_R20 );
  PLA_Obj_local_length       ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 && value_R12 != PLA_NO_DIMENSION ) {
    PLA_Warning( "inconsistent local_length R20 vs R12" );
    printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_local_width_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_LOCAL_WIDTH      ( obj, &value_R20 );
  PLA_Obj_local_width       ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 && value_R12 != PLA_NO_DIMENSION ) {
    PLA_Warning( "inconsistent local_width R20 vs R12" );
    printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_local_ldim_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_LOCAL_LDIM         ( obj, &value_R20 );
  PLA_Obj_local_ldim          ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 ) {
    if ( value_R20 != 0 || value_R12 != PLA_NO_DIMENSION ){
      PLA_Warning( "inconsistent local_ldim R20 vs R12" );
      printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
      return_value = FALSE;
    }
  }

  return return_value;
}

int PLA_Obj_local_stride_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    value_R20, value_R12,
    me;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_LOCAL_STRIDE       ( obj, &value_R20 );
  PLA_Obj_local_stride        ( obj->pla_obj, &value_R12 );
  if ( value_R20 != value_R12 ) {
    PLA_Warning( "inconsistent local_stride R20 vs R12" );
    printf("%3d: %d vs %d\n", me, value_R20, value_R12 );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_template_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    me;

  PLA_Template 
    templ_R20, templ_R12;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_TEMPLATE           ( obj, &templ_R20 );
  PLA_Obj_template            ( obj->pla_obj, &templ_R12 );
  if ( templ_R20 != templ_R12 ) {
    PLA_Warning( "inconsistent template R20 vs R12" );
    return_value = FALSE;
  }

  return return_value;
}

int PLA_Obj_sizes_and_owners_consistent_with_R12( PLA_Obj obj )
{
  int 
    return_value = TRUE,
    me,
    size_R20, size_R12,
    owner_R20, owner_R12;

  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  PLA_OBJ_SIZE_TOP           ( obj, &size_R20 );
  PLA_OBJ_OWNER_TOP          ( obj, &owner_R20 );
  PLA_Obj_split_size ( obj->pla_obj, PLA_SIDE_TOP, &size_R12, &owner_R12 );
  if ( size_R20 != size_R12 ) {
    PLA_Warning( "inconsistent size_top R20 vs R12" );
    printf("%3d: %d vs %d\n", me, size_R20, size_R12 );
    return_value = FALSE;
  }

  if ( owner_R20 != owner_R12 ) {
    if ( !( 
	   ( owner_R20 == PLA_ALL_ROWS && owner_R12 == PLA_ALL_ROWS ) ||
	   ( owner_R20 == PLA_ALL_COLS && owner_R12 == PLA_ALL_COLS ) ) ) {
      PLA_Warning( "inconsistent owner_top R20 vs R12" );
      printf("%3d: %d vs %d\n", me, owner_R20, owner_R12 );
      return_value = FALSE;
    }
  }

  PLA_OBJ_SIZE_BOTTOM        ( obj, &size_R20 );
  PLA_OBJ_OWNER_BOTTOM       ( obj, &owner_R20 );
  PLA_Obj_split_size ( obj->pla_obj, PLA_SIDE_BOTTOM, &size_R12, &owner_R12 );
  if ( size_R20 != size_R12 ) {
    PLA_Warning( "inconsistent size_bottom R20 vs R12" );
    printf("%3d: %d vs %d\n", me, size_R20, size_R12 );
    return_value = FALSE;
  }
  if ( owner_R20 != owner_R12 ) {
    if ( !( 
	   ( owner_R20 == PLA_ALL_ROWS && owner_R12 == PLA_ALL_ROWS ) ||
	   ( owner_R20 == PLA_ALL_COLS && owner_R12 == PLA_ALL_COLS ) ) ) {
      int length; 

      PLA_OBJ_GLOBAL_LENGTH( obj, &length );
      if ( length != 0 ) {
	PLA_Warning( "inconsistent owner_bottom R20 vs R12" );
	printf("%3d: %d vs %d\n", me, owner_R20, owner_R12 );
	return_value = FALSE;
      }
    }
  }

  PLA_OBJ_SIZE_LEFT           ( obj, &size_R20 );
  PLA_OBJ_OWNER_LEFT          ( obj, &owner_R20 );
  PLA_Obj_split_size ( obj->pla_obj, PLA_SIDE_LEFT, &size_R12, &owner_R12 );
  if ( size_R20 != size_R12 ) {
    PLA_Warning( "inconsistent size_left R20 vs R12" );
    printf("%3d: %d vs %d\n", me, size_R20, size_R12 );
    return_value = FALSE;
  }
  if ( owner_R20 != owner_R12 ) {
    if ( !( 
	   ( owner_R20 == PLA_ALL_ROWS && owner_R12 == PLA_ALL_ROWS ) ||
	   ( owner_R20 == PLA_ALL_COLS && owner_R12 == PLA_ALL_COLS ) ) ) {
      PLA_Warning( "inconsistent owner_left R20 vs R12" );
      printf("%3d: %d vs %d\n", me, owner_R20, owner_R12 );
      return_value = FALSE;
    }
  }

  PLA_OBJ_SIZE_RIGHT           ( obj, &size_R20 );
  PLA_OBJ_OWNER_RIGHT          ( obj, &owner_R20 );
  PLA_Obj_split_size ( obj->pla_obj, PLA_SIDE_RIGHT, &size_R12, &owner_R12 );
  if ( size_R20 != size_R12 ) {
    PLA_Warning( "inconsistent size_right R20 vs R12" );
    printf("%3d: %d vs %d\n", me, size_R20, size_R12 );
    return_value = FALSE;
  }
  if ( owner_R20 != owner_R12 ) {
    if ( !( 
	   ( owner_R20 == PLA_ALL_ROWS && owner_R12 == PLA_ALL_ROWS ) ||
	   ( owner_R20 == PLA_ALL_COLS && owner_R12 == PLA_ALL_COLS ) ) ) {
      int width;

      PLA_OBJ_GLOBAL_WIDTH( obj, &width );
      if ( width != 0 ) {
	PLA_Warning( "inconsistent owner_right R20 vs R12" );
	printf("%3d: %d vs %d\n", me, owner_R20, owner_R12 );
	return_value = FALSE;
      }
    }
  }

  return return_value;
}
#endif

