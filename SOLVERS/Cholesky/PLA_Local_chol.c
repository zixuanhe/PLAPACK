/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

double sqrt( double );

int PLA_Local_chol ( int uplo, PLA_Obj A )
{
  int 
    value = PLA_SUCCESS;

  /* Perform parameter and error checking */
  /*  if ( PLA_ERROR_CHECKING )    
      value = PLA_Local_chol_enter( uplo, A ); */

  if ( value == 0 ){
    switch ( uplo ){
    case PLA_LOWER_TRIANGULAR: 
      value = PLA_Local_chol_lower( A );
      break;
    case PLA_UPPER_TRIANGULAR: 
      /*    value = PLA_Local_chol_upper( A ); */
      value = -1;
    }
  }

  /* Perform parameter and error checking */
  /*  if ( PLA_ERROR_CHECKING )    
      value = PLA_Local_chol_exit( uplo, A ); */

  return( value );
}
  

int PLA_Local_chol_lower ( PLA_Obj A )
{
   int 
     value,
     local_length, local_width, local_ldim;

   void 
     *buffer;

   MPI_Datatype 
     datatype;

   PLA_Obj_local_length( A, &local_length );
   PLA_Obj_local_width( A, &local_width );
   if( local_length == 0 || local_width == 0 ) 
     return PLA_SUCCESS;

   PLA_Obj_local_ldim  ( A, &local_ldim );
   PLA_Obj_datatype    ( A, &datatype );
   PLA_Obj_local_buffer( A, &buffer );

   if ( datatype == MPI_DOUBLE )
     value = PLA_dchol ( local_length, ( double * ) buffer, local_ldim );
   /*
     PLA_dpotrf( "L", &local_length, buffer, &local_ldim, &value );
   */
   else if ( datatype == MPI_FLOAT )
     value = PLA_schol ( local_length, ( float * ) buffer,  local_ldim ); 
   /*
     PLA_spotrf( "L", &local_length, buffer, &local_ldim, &value );
   */
   else if ( datatype == MPI_DOUBLE_COMPLEX )
     value = PLA_zchol ( local_length, ( float * ) buffer,  local_ldim );
     /*
     PLA_zpotrf( "L", &local_length, buffer, &local_ldim, &value );
     */
   else if ( datatype == MPI_COMPLEX )
     value = PLA_cchol ( local_length, ( float * ) buffer,  local_ldim );
   /*
     PLA_cpotrf( "L", &local_length, buffer, &local_ldim, &value );
   */
   else
     value = PLA_Warning( "Invalid datatype in PLA_Local_chol" );

   return value;
}


/************************************************
  Note on error return codes.

  Return code of 0 means normal completion.

  Return code of i (i > 0), means non-positive
  diagonal element was encountered at the i_th
  diagonal entry (the first entry is '1' rather 
  than '0' to distinguish from normal return.
***********************************************/


int PLA_dchol ( int n, double * buffer, int lld )
{
   int k,
       length,
       int_one = 1;
   double invert,
          min_one = -1.0e+0,
          one = 1.0e+0;

#define a(i,j) buffer [ ( j ) * lld + ( i ) ]

   for ( k = 0 ; k < n ; k++ )
   {
      if ( k > 0 )
      {
         length = n - k;
         PLA_dgemv ( "N", & length, & k,
                     & min_one, & a ( k, 0 ), & lld, & a ( k, 0 ), & lld,
                     & one, & a ( k, k ), & int_one );
      }

      if ( a (k,k) <= 0.0 ) {
	printf("a(k,k) = %lf\n", a(k,k) );
	return (k+1);
      }

      a ( k, k ) = sqrt ( a ( k, k ) );

      length = n - k - 1;

      if ( length > 0 )
      {
         invert = one / a ( k, k );
         PLA_dscal ( & length, & invert, & a ( k + 1, k ), & int_one );
      }
   }

#undef a

   return ( 0 );
}

int PLA_schol ( int n, float * buffer, int lld )
{
   int k,
       length,
       int_one = 1;
   float invert,
         min_one = -1.0e+0,
         one = 1.0e+0;

#define a(i,j) buffer [ ( j ) * lld + ( i ) ]

   for ( k = 0 ; k < n ; k++ )
   {

      if ( k > 0 )
      {
         length = n - k;
         PLA_sgemv ( "N", & length, & k,
                     & min_one, & a ( k, 0 ), & lld, & a ( k, 0 ), & lld,
                     & one, & a ( k, k ), & int_one );
      }

      if ( a ( k, k) <= 0.0 ) return (k+1);

      a ( k, k ) = (float) sqrt ( (double) a ( k, k ) );

      length = n - k - 1;

      if ( length > 0 )
      {
         invert = one / a ( k, k );
         PLA_sscal ( & length, & invert, & a ( k + 1, k ), & int_one );
      }
   }

#undef a

   return ( 0 );
}


int PLA_cchol ( int n, PLA_COMPLEX * buffer, int lld )
{
  int j, i_one = 1, length;

  PLA_COMPLEX 
    invert = {0.0, 0.0},
    min_one = {-1.0,0.0};

#define a(i,j) buffer [ ( j ) * lld + ( i ) ]

  for ( j=0; j<n; j++ ){
    a(j,j).real = (float) sqrt( (double) a(j,j).real );
    a(j,j).imaginary = 0.0;
    length = n-j-1;
    invert.real = 1.0/a(j,j).real;
    PLA_cscal( &length, &invert, &a(j+1,j), &i_one );
    PLA_cher( "L", &length, &min_one, &a(j+1,j  ), &i_one, 
                                         &a(j+1,j+1), &lld );
   }

  return PLA_SUCCESS;
#undef a
}

int PLA_zchol ( int n, PLA_DOUBLE_COMPLEX * buffer, int lld )
{
  int j, i_one = 1, length;

  PLA_DOUBLE_COMPLEX 
    invert,
    min_one = {-1.0, 0}, one = {1.0, 0};

#define a(i,j) buffer [ ( j ) * lld + ( i ) ]

  invert.imaginary = 0.0;

  for ( j=0; j<n; j++ ){
    a(j,j).real = sqrt( a(j,j).real );
    a(j,j).imaginary = 0.0;
    length = n-j-1;
    invert.real = 1.0/a(j,j).real;
    PLA_zscal( &length, &invert, &a(j+1,j), &i_one );
    PLA_zher( "L", &length, &min_one, &a(j+1,j  ), &i_one, 
                                       &a(j+1,j+1), &lld );
   }

#undef a

  return( PLA_SUCCESS );
}

     


