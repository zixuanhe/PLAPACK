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

/*
#if MANUFACTURE == CRAY
#define lseek      ffseek
#define reada      ffreada
#define writea     ffwritea
#define listio     fflistio
#define listreq    fflistreq
#define iosw       ffsw
#endif
*/

int PLA_Read_async( int size, int fd, long offset, void *buffer,
#if MANUFACTURE == CRAY
		     struct ffsw *status)
#else
                     int *status )
#endif
{
  int 
    i, temp;

#if MANUFACTURE == CRAY
  if ( ( temp = ffseek( fd, offset, 0 ) ) != offset ){
    PLA_Abort( "error seeking in file before reading", __LINE__, __FILE__ );
  }

  /* call ffreada */

  if ( ( temp = ffreada( fd, buffer, size, status ) ) <0 ){
    printf("size=%d temp=%d errno=%d\n", size, temp, status->sw_error);
    PLA_Warning( "error reading file" );
  }
#else
  PLA_Abort( "function not yet implemented", __LINE__, __FILE__ );
#endif

  return PLA_SUCCESS;
}

int PLA_Write_async( int size, void *buffer, int fd, long offset,
#if MANUFACTURE == CRAY
		      struct ffsw *status)
#else
                      int *status )
#endif
{
  int 
    i;


#if MANUFACTURE == CRAY
  ffseek( fd, offset, 0 );

  /* call writea */

  if ( ffwritea( fd, buffer, size, status )  <0)
    PLA_Abort( "error writing file", __LINE__, __FILE__ );
#else
  PLA_Abort( "function not yet implemented", __LINE__, __FILE__ );
#endif

  return PLA_SUCCESS;
}


int PLA_Read_matrix_async( int length, int width, 
		      int fd, int buffer_from, int ldim_from,
		              char *buffer_to,   int ldim_to,
#if MANUFACTURE == CRAY
		      struct ffsw *status)
#else
                      int *status )
#endif
{
  int j;
#if MANUFACTURE == CRAY
  struct fflistreq request;
  

  if ( ldim_from == length && ldim_to == length )
    PLA_Read_async( ( long ) ( length*width ), fd, buffer_from, 
	           buffer_to, status );
  else{
    printf("LISTIO READ\n");
    /* Call listio */
    request.li_opcode = LO_READ;
    request.li_flags = LF_LSEEK;
    request.li_offset = buffer_from;
    request.li_fildes = fd;
    request.li_buf = buffer_to;
    request.li_nbyte = length;
    request.li_status = status;
    request.li_signo = 0;
    request.li_nstride = width;
    request.li_memstride = ldim_to;
    request.li_filstride = ldim_from;

    fflistio(LC_START, &request, 1);
  } 
#else
  PLA_Abort( "function not yet implemented", __LINE__, __FILE__ );
#endif
}


int PLA_Write_matrix_async( int length, int width, 
		       char *buffer_from,   int ldim_from,
		       int fd, int buffer_to, int ldim_to,
#if MANUFACTURE == CRAY
		       struct ffsw *status)
#else
                      int *status )
#endif
{
  int j;
#if MANUFACTURE == CRAY
  struct fflistreq request;

  if ( ldim_from == length && ldim_to == length )
    PLA_Write_async( ( long ) ( length*width ), buffer_from, fd,
	           buffer_to, status );
  else{
    printf("LISTIO WRITE\n");
    /* Call listio */
    
    request.li_opcode = LO_WRITE;
    request.li_flags = LF_LSEEK;
    request.li_offset = buffer_to;
    request.li_fildes = fd;
    request.li_buf = buffer_from;
    request.li_nbyte = length;
    request.li_status = status;
    request.li_signo = 0;
    request.li_nstride = width;
    request.li_memstride = ldim_from;
    request.li_filstride = ldim_to;

    fflistio(LC_START, &request, 1);

  } 
#else
  PLA_Abort( "function not yet implemented", __LINE__, __FILE__ );
#endif
}




