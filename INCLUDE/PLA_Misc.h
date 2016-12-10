/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

typedef struct PLA_COMPLEX{
  float real, imaginary;
} PLA_COMPLEX;

typedef struct PLA_DOUBLE_COMPLEX{
  double real, imaginary;
} PLA_DOUBLE_COMPLEX;

#ifndef min
#define min(x,y) ( (x) < (y) ? (x) : (y) )
#endif

#ifndef dabs
#define dabs(x) ( (x) < 0 ? -(x) : (x) )
#endif

#ifndef MPI_COMPLEX
  #if MANUFACTURE != SGI && ! (MANUFACTURE == CRAY && MACHINE_TYPE == CRAYPVP)
    MPI_Datatype MPI_COMPLEX;
    #define PLA_MPI_COMPLEX TRUE
  #endif
#endif

#ifndef MPI_DOUBLE_COMPLEX
  #if MANUFACTURE != SGI && ! (MANUFACTURE == CRAY && MACHINE_TYPE == CRAYPVP)
    MPI_Datatype MPI_DOUBLE_COMPLEX;
    #define PLA_MPI_DOUBLE_COMPLEX TRUE
  #endif
#endif

double drand48();

#define PLA_METHOD_INV 0
#define PLA_METHOD_FACTORS 1
#define PLA_METHOD_STABLE 2
#define PLA_METHOD_TRSM 3

