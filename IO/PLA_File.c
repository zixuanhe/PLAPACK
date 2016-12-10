/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include <sys/types.h>
#include <sys/stat.h>

#if MANUFACTURE == CRAY

#define open       ffopen
#define close      ffclose
#define lseek      ffseek
#define read       ffread
#define write      ffwrite

#endif

int PLA_File_create( char *filename, MPI_Datatype datatype, int size, int *fd )
{
  char
    message[40];

  int 
    typesize;


  if (( *fd = open( filename, O_RDWR|O_CREAT|O_TRUNC, 0700 )) == -1 ){
    sprintf( message, "Cannot open file %s", filename );
    PLA_Abort( message, __LINE__, __FILE__ );
  }

  MPI_Type_size( datatype, &typesize );

  lseek( *fd, ( long ) ( size * typesize ), 0 );

  if ( write( *fd, "\0", 1) != size )
    
  return PLA_SUCCESS;
}

int PLA_File_open( char *filename, MPI_Datatype datatype, int size, int *fd )
{
  char
    message[40];
 
  if ( -1 == ( *fd = open( filename, O_RDWR ) ) ){
    sprintf( message, "Cannot open file %s", filename );
    PLA_Abort( message, __LINE__, __FILE__ );
  }

  return PLA_SUCCESS;
}

int PLA_File_close( int fd )
{
  char
    message[40];

  if ( -1 == close( fd ) ){
    sprintf( message, "Trouble closing file" );
    PLA_Abort( message, __LINE__, __FILE__ );
  }

  return PLA_SUCCESS;
}

int PLA_File_unlink( char *filename )
{
  char
    message[40];

  if ( -1 == unlink( filename ) ){
    sprintf( message, "Trouble unlinking file %s", filename );
    PLA_Abort( message, __LINE__, __FILE__ );
  }

  return PLA_SUCCESS;
}

int PLA_Read( int size, int fd, long offset, void *buffer )
{
  int 
    i, temp;

  double
    time;

  time = MPI_Wtime();
  if ( ( temp = lseek( fd, offset, 0 ) ) != offset ){
    PLA_Abort( "error seeking in file before reading", __LINE__, __FILE__ );
  }

  if ( ( temp = read( fd, buffer, size ) ) != size ){
    PLA_Warning( "error reading file" );
  }
  time = MPI_Wtime() - time;

  PLA_File_stats_read(size, time);

  return PLA_SUCCESS;
}

int PLA_Write( int size, void *buffer, int fd, long offset )
{
  int 
    i;

  double
    time;

  time = MPI_Wtime();
  lseek( fd, offset, 0 );

  if ( write( fd, buffer, size ) != size )
    PLA_Abort( "error writing file", __LINE__, __FILE__ );
  time = MPI_Wtime() - time;
  PLA_File_stats_write(size, time);
  return PLA_SUCCESS;
}


int PLA_Read_matrix( int length, int width, 
		      int fd, int buffer_from, int ldim_from,
		              char *buffer_to,   int ldim_to )
{
  int j;

  if ( ldim_from == length && ldim_to == length )
    PLA_Read( ( long ) ( length*width ), fd, buffer_from, 
	           buffer_to );
  else{

    for ( j=0; j<width; j++ ){
      PLA_Read( ( long ) ( length ), fd, buffer_from,
		     buffer_to );
      buffer_from += ldim_from;
      buffer_to +=   ldim_to;
    } 

  } 
}


int PLA_Write_matrix( int length, int width, 
		       char *buffer_from,   int ldim_from,
		       int fd, int buffer_to, int ldim_to )
{
  int j;

  if ( ldim_from == length && ldim_to == length )
    PLA_Write( ( long ) ( length*width ), buffer_from, fd,
	           buffer_to );
  else{

    for ( j=0; j<width; j++ ){
      PLA_Write( ( long ) ( length ), buffer_from, fd, 
		     buffer_to );
      buffer_from += ldim_from;
      buffer_to +=   ldim_to;
    } 

  } 
}
