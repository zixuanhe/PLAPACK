/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

#if MANUFACTURE == IBM
#define MPI_Fint int
#define MPI_Comm_c2f(obj) (obj)
#define MPI_Comm_f2c(obj) (obj)
#define MPI_Op_f2c(op)    (op)
#endif


#if MANUFACTURE == CRAY
#define MPI_Fint int
#define MPI_Comm_c2f(obj) (obj)
#define MPI_Comm_f2c(obj) (obj)
#define MPI_Op_f2c(op)    (op)
#endif

#if MANUFACTURE == PC
#undef PLA_FOR2C
#define PLA_FOR2C( name ) name ## __ 
#endif

#define NUM_OF_TYPES 6

static int initialized = FALSE;
static MPI_Datatype types_C[ NUM_OF_TYPES ];
static int types_F[ NUM_OF_TYPES ];

pla_type_conversion_init( MPI_Datatype *type_list )
{
  type_list[0] = MPI_CHAR;
  type_list[1] = MPI_INT;
  type_list[2] = MPI_FLOAT;
  type_list[3] = MPI_DOUBLE;
  type_list[4] = MPI_COMPLEX;
  type_list[5] = MPI_DOUBLE_COMPLEX;
}

int pla_type_c2f( MPI_Datatype datatype )
{
  int i;

  if (!initialized) {
    pla_type_conversion_init( types_C );
#if MANUFACTURE == CRAY
    PLA_TYPE_CONVERSION_CRAY_INIT_F( types_F );
#else
    PLA_FOR2C( pla_type_conversion_init_f )( types_F );
#endif    
    initialized = TRUE;
  }

  
  for ( i=0; i<NUM_OF_TYPES; i++ ) 
    if ( datatype == types_C[ i ] ) break;
    
  if ( i >= NUM_OF_TYPES ) {
    printf("error in PLA_Type_c2f: type not found\n");
    exit( 0 );
  }

  return types_F[i];
}

MPI_Datatype pla_type_f2c( int datatype )
{
  int i;

  if (!initialized) {
    pla_type_conversion_init( types_C );
#if MANUFACTURE == CRAY
    PLA_TYPE_CONVERSION_CRAY_INIT_F( types_F );
#else
    PLA_FOR2C( pla_type_conversion_init_f )( types_F );
#endif    
    
    initialized = TRUE;
  }

  
  for ( i=0; i<NUM_OF_TYPES; i++ ) 
    if ( datatype == types_F[ i ] ) break;
    
  if ( i >= NUM_OF_TYPES ) {
    printf("error in PLA_Type_f2c: type not found\n");
    exit( 0 );
  }
  
  return types_C[ i ];
}


