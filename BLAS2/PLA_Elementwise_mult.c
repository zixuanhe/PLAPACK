/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

#define A(i,j) local_buff_A[ (j)*(local_ldim_A) + (i) ]
#define B(i,j) local_buff_B[ (j)*(local_ldim_B) + (i) ]
#define C(i,j) local_buff_C[ (j)*(local_ldim_C) + (i) ]

int PLA_Elementwise_mult( PLA_Obj A, PLA_Obj B, PLA_Obj C )
{
  int 
    objtype_A, objtype_B, objtype_C,
    local_length_A, local_width_A,
    local_length_B, local_width_B,
    local_length_C, local_width_C,
    global_length_A, global_width_A,
    global_length_B, global_width_B,
    global_length_C, global_width_C,
    global_align_row_A, global_align_col_A,
    global_align_row_B, global_align_col_B,
    global_align_row_C, global_align_col_C,
    local_ldim_A, local_ldim_B, local_ldim_C,
    i, j; 

  MPI_Datatype
    datatype;
  
  PLA_Obj_objtype( A, &objtype_A );
  PLA_Obj_objtype( B, &objtype_B );
  PLA_Obj_objtype( C, &objtype_C );

  if ( objtype_A != objtype_B || objtype_A != objtype_C )
    PLA_Abort( "PLA_Elementwise_mult: objtype must all match", __LINE__, __FILE__ );

  PLA_Obj_global_length( A, &global_length_A );
  PLA_Obj_global_width ( A, &global_width_A );
  PLA_Obj_global_length( B, &global_length_B );
  PLA_Obj_global_width ( B, &global_width_B );
  PLA_Obj_global_length( C, &global_length_C );
  PLA_Obj_global_width ( C, &global_width_C );

  if ( global_length_A != global_length_B || global_length_A != global_length_C )
    PLA_Abort( "PLA_Elementwise_mult: row dimensions must all match", __LINE__, __FILE__ );

  if ( global_width_A != global_width_B || global_width_A != global_width_C )
    PLA_Abort( "PLA_Elementwise_mult: column dimensions must all match", __LINE__, __FILE__ );

  if ( objtype_A == PLA_MATRIX ){
    PLA_Obj_global_align_row( A, &global_align_row_A );
    PLA_Obj_global_align_row( B, &global_align_row_B );
    PLA_Obj_global_align_row( C, &global_align_row_C );
  
    if ( global_align_row_A != global_align_row_B || global_align_row_A != global_align_row_C )
      PLA_Abort( "PLA_Elementwise_mult: row alignments must all match", __LINE__, __FILE__ );

    PLA_Obj_global_align_col( A, &global_align_col_A );
    PLA_Obj_global_align_col( B, &global_align_col_B );
    PLA_Obj_global_align_col( C, &global_align_col_C );
  
    if ( global_align_col_A != global_align_col_B || global_align_col_A != global_align_col_C )
      PLA_Abort( "PLA_Elementwise_mult: col alignments must all match", __LINE__, __FILE__ );
  }
  else if ( objtype_A == PLA_MVECTOR ){
    PLA_Obj_global_align_row( A, &global_align_row_A );
    PLA_Obj_global_align_row( B, &global_align_row_B );
    PLA_Obj_global_align_row( C, &global_align_row_C );
  
    if ( global_align_row_A != global_align_row_B || global_align_row_A != global_align_row_C )
      PLA_Abort( "PLA_Elementwise_mult: row alignments must all match", __LINE__, __FILE__ );
  }
  else 
      PLA_Abort( "PLA_Elementwise_mult: objtype not yet supported ", __LINE__, __FILE__ );

  PLA_Obj_local_length( A, &local_length_A );
  PLA_Obj_local_width ( A, &local_width_A );
    
  PLA_Obj_local_ldim( A, &local_ldim_A );
  PLA_Obj_local_ldim( B, &local_ldim_B );
  PLA_Obj_local_ldim( C, &local_ldim_C );

  PLA_Obj_datatype( A, &datatype );
  
  if ( MPI_DOUBLE == datatype ){
    double 
      *local_buff_A, *local_buff_B, *local_buff_C;
    
    PLA_Obj_local_buffer( A, ( void **) &local_buff_A );
    PLA_Obj_local_buffer( B, ( void **) &local_buff_B );
    PLA_Obj_local_buffer( C, ( void **) &local_buff_C );
  
    for ( j=0; j<local_width_A; j++ )
      for ( i=0; i<local_length_A; i++ )

	C( i, j ) = A( i, j ) * B( i, j );
  }
  else if ( MPI_FLOAT == datatype ){
    float 
      *local_buff_A, *local_buff_B, *local_buff_C;
    
    PLA_Obj_local_buffer( A, ( void **) &local_buff_A );
    PLA_Obj_local_buffer( B, ( void **) &local_buff_B );
    PLA_Obj_local_buffer( C, ( void **) &local_buff_C );
  
    for ( j=0; j<local_width_A; j++ )
      for ( i=0; i<local_length_A; i++ )

	C( i, j ) = A( i, j ) * B( i, j );
  }
  else 
    PLA_Abort( "PLA_Elementwise_mult: datatype not yet supported ", __LINE__, __FILE__ );
}
