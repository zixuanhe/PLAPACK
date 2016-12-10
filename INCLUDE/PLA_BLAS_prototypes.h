/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/*
 * Single Precision
 */

void  srotg ( float *, float *, float *, float * ); 
void  srot  ( int *, float *, int *, float *, int *, 
	                             float *, float * );
void  sswap ( int *, float *, int *, float *, int * );
void  scopy ( int *, float *, int *, float *, int * );
void  saxpy ( int *, float *, float *, int *, float *, int * );
void  sscal ( int *, float *, float *, int * );
float sdot ( int *, float *, int *, float *, int * );
float snrm2 ( int *, float *, int * );
float sasum ( int *, float *, int * );
int   isamax ( int *, float *, int * );
void  slacpy ( char *, int *, int *, float *, int *, float *, int * );

/*
 * Double Precision
 */

void  drotg ( double *, double *, double *, double * ); 
void  drot  ( int *, double *, int *, double *, int *, 
	                             double *, double * );
void  dswap ( int *, double *, int *, double *, int * );
void  dcopy ( int *, double *, int *, double *, int * );
void  daxpy ( int *, double *, double *, int *, double *, int * );
void  dscal ( int *, double *, double *, int * );
double ddot ( int *, double *, int *, double *, int * );
double dnrm2 ( int *, double *, int * );
double dasum ( int *, double *, int * );
int   idamax ( int *, double *, int * );
void  dlacpy ( char *, int *, int *, double *, int *, double *, int * );

/*
 * Single Precision Complex
 */

void  crotg ( PLA_COMPLEX *, PLA_COMPLEX *, PLA_COMPLEX *, PLA_COMPLEX * ); 
void  crot  ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int *, 
	                             PLA_COMPLEX *, PLA_COMPLEX * );
void  cswap ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
void  ccopy ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
void  caxpy ( int *, PLA_COMPLEX *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
void  cscal ( int *, PLA_COMPLEX *, PLA_COMPLEX *, int * );
PLA_COMPLEX cdot ( int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );
float cnrm2 ( int *, PLA_COMPLEX *, int * );
float casum ( int *, PLA_COMPLEX *, int * );
int   icamax ( int *, PLA_COMPLEX *, int * );
void  clacpy ( char *, int *, int *, PLA_COMPLEX *, int *, PLA_COMPLEX *, int * );

/*
 * Double Precision Complex
 */

void  zrotg ( PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX * ); 
void  zrot  ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int *, 
	                             PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX * );
void  zswap ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
void  zcopy ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
void  zaxpy ( int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
void  zscal ( int *, PLA_DOUBLE_COMPLEX *, PLA_DOUBLE_COMPLEX *, int * );
PLA_DOUBLE_COMPLEX zdot ( int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );
double znrm2 ( int *, PLA_DOUBLE_COMPLEX *, int * );
double zasum ( int *, PLA_DOUBLE_COMPLEX *, int * );
int   izamax ( int *, PLA_DOUBLE_COMPLEX *, int * );
void  zlacpy ( char *, int *, int *, PLA_DOUBLE_COMPLEX *, int *, PLA_DOUBLE_COMPLEX *, int * );

