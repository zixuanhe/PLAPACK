/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

static int PLA_DUMMY_ADRESS = 1;

/* #if MANUFACTURE == CRAY */
  #define PLA_CRAY_NULL_FIX( name ) ( ( name ) ? ( name ) : 1 )
  #define BF(buffer) (NULL==buffer) ? &PLA_DUMMY_ADRESS : buffer

/*
  #define PLA_MPI_Allgatherv(a,b,c,d,e,f,g,h)   \
       do {                                     \
          MPI_Barrier( h );                     \ 
	  MPI_Allgatherv(a,b,c,d,e,f,g,h);      \
       } while(0)
*/
  #define PLA_MPI_Allgatherv(a,b,c,d,e,f,g,h)   \
       do {                                     \
	  MPI_Allgatherv(a,b,c,d,e,f,g,h);      \
       } while(0)

/* #else
  #define PLA_CRAY_NULL_FIX( name )  name
  #define BF(buffer) buffer 
  #define PLA_MPI_Allgatherv(a,b,c,d,e,f,g,h)   \
	  MPI_Allgatherv    (a,b,c,d,e,f,g,h) 
#endif */

