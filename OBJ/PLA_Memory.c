/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/***************************************************************************/

static total_size_malloced = 0;
static max_size_malloced = 0;

void * PLA_malloc ( size_t size )

/*----------------------------------------------------------------------------

Purpose : Create allocated space


IN     size              size of buffer to be created

----------------------------------------------------------------------------*/
{
  void *buffer;

  if ( size < 0 ) {
    printf( "attempt to malloc negative number of bytes\n");
    exit( 0 );
  }

  buffer = ( void * )malloc( size );

  if ( buffer == NULL && size > 0 ){
    printf( "malloc failed in PLA_malloc \n");
    exit( 0 );
  }

  total_size_malloced += malloc_howbig( buffer ); 

  if ( total_size_malloced > max_size_malloced )
    max_size_malloced = total_size_malloced;

  return buffer;
}

void *PLA_calloc( size_t size, size_t typesize )
{
  void *return_value;
  char *tempp;
  int i;

  return_value = PLA_malloc( (size_t) size*typesize );
  tempp = ( char * ) return_value;

  for ( i=0; i<size*typesize; i++ ) *tempp++ = 0;

  return return_value;
}


/***************************************************************************/

void * PLA_free ( void *buffer )

/*----------------------------------------------------------------------------

Purpose : Free allocated space


IN     buffer           address of buffer to be freed

----------------------------------------------------------------------------*/
{
  if ( buffer != NULL ) {
    total_size_malloced -= malloc_howbig( buffer );
    free( buffer );
  }
}

/***************************************************************************/

int PLA_Total_size_malloced ( )

/*----------------------------------------------------------------------------

Purpose : Return total size of malloced space

----------------------------------------------------------------------------*/
{
  return total_size_malloced;
}

/***************************************************************************/

int PLA_Max_size_malloced ( )

/*----------------------------------------------------------------------------

Purpose : Return max size of malloced space to date

----------------------------------------------------------------------------*/
{
  return max_size_malloced;
}


#if MANUFACTURE == CRAY
#else  
int malloc_howbig( void *tempp )
{
  return 0;
}
#endif









