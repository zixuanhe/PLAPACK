/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

#define A(i,j,lda) a[ ((j)-1)*(lda) + (i)-1 ]
#define Piv(i) pivot[i-1]

PLA_DOUBLE_COMPLEX double_complex_inverse( PLA_DOUBLE_COMPLEX );
PLA_COMPLEX complex_inverse( PLA_COMPLEX );

/* double precision */

void PLA_dlu( int *m, int *n, double *a, int *lda, int *pivot )
{
  int
    j, index, m_min_j_p_1, m_min_j, n_min_j, i_one = 1, i_four = 4;
  double
    ajj_inv, minus_one = -1.0;

  for ( j=1; j<=min( *n, *m ); j++ ){ 
    m_min_j_p_1 = *m-j+1;
    m_min_j = *m-j;
    n_min_j = *n-j;

    index = PLA_idamax( &m_min_j_p_1, &A(j,j,*lda), &i_one );

    Piv(j) = index-1;

    PLA_dswap( n, &A(j,1,*lda), lda, &A(j+index-1,1,*lda), lda ); 

    ajj_inv = 1.0/A(j,j,*lda);

    PLA_dscal( &m_min_j, &ajj_inv, &A(j+1,j,*lda), &i_one );

    PLA_dger( &m_min_j, &n_min_j, &minus_one, &A(j+1,j,*lda), &i_one,
	       &A(j,j+1,*lda), lda, &A(j+1,j+1,*lda ), lda );
  }
}

/* single precision */

void PLA_slu( int *m, int *n, float *a, int *lda, int *pivot )
{
  int
    j, index, m_min_j_p_1, m_min_j, n_min_j, i_one = 1, i_four = 4;
  float
    ajj_inv, minus_one = -1.0;

  for ( j=1; j<=min( *n, *m ); j++ ){ 
    m_min_j_p_1 = *m-j+1;
    m_min_j = *m-j;
    n_min_j = *n-j;

    index = PLA_isamax( &m_min_j_p_1, &A(j,j,*lda), &i_one );

    Piv(j) = index-1;

    PLA_sswap( n, &A(j,1,*lda), lda, &A(j+index-1,1,*lda), lda ); 

    ajj_inv = 1.0/A(j,j,*lda);

    PLA_sscal( &m_min_j, &ajj_inv, &A(j+1,j,*lda), &i_one );
    
    PLA_sger( &m_min_j, &n_min_j, &minus_one, &A(j+1,j,*lda), &i_one,
	       &A(j,j+1,*lda), lda, &A(j+1,j+1,*lda ), lda );
  }
}

/* double precision complex */

void PLA_zlu( int *m, int *n, PLA_DOUBLE_COMPLEX *a, int *lda, int *pivot )
{
  int
    j, index, m_min_j_p_1, m_min_j, n_min_j, i_one = 1, i_four = 4;
  PLA_DOUBLE_COMPLEX
    ajj_inv, minus_one = {-1.0, 0.0};

  for ( j=1; j<=min( *n, *m ); j++ ){ 
    m_min_j_p_1 = *m-j+1;
    m_min_j = *m-j;
    n_min_j = *n-j;

    index = PLA_izamax( &m_min_j_p_1, &A(j,j,*lda), &i_one );

    Piv(j) = index-1;

    PLA_zswap( n, &A(j,1,*lda), lda, &A(j+index-1,1,*lda), lda ); 

    ajj_inv = double_complex_inverse( A(j,j,*lda) );

    PLA_zscal( &m_min_j, &ajj_inv, &A(j+1,j,*lda), &i_one );
    
    PLA_zger( &m_min_j, &n_min_j, &minus_one, &A(j+1,j,*lda), &i_one,
	       &A(j,j+1,*lda), lda, &A(j+1,j+1,*lda ), lda );
  }
}

/* precision complex */

void PLA_clu( int *m, int *n, PLA_COMPLEX *a, int *lda, int *pivot )
{
  int
    j, index, m_min_j_p_1, m_min_j, n_min_j, i_one = 1, i_four = 4;
  PLA_COMPLEX
    ajj_inv, minus_one = {-1.0, 0.0};

  for ( j=1; j<=min( *n, *m ); j++ ){ 
    m_min_j_p_1 = *m-j+1;
    m_min_j = *m-j;
    n_min_j = *n-j;

    index = PLA_icamax( &m_min_j_p_1, &A(j,j,*lda), &i_one );

    Piv(j) = index-1;

    PLA_cswap( n, &A(j,1,*lda), lda, &A(j+index-1,1,*lda), lda ); 

    ajj_inv = complex_inverse( A(j,j,*lda) );

    PLA_cscal( &m_min_j, &ajj_inv, &A(j+1,j,*lda), &i_one );
    
    PLA_cger( &m_min_j, &n_min_j, &minus_one, &A(j+1,j,*lda), &i_one,
	       &A(j,j+1,*lda), lda, &A(j+1,j+1,*lda ), lda );
  }
}


#define abs_max(x,y) ( dabs(x) > dabs(y) ? dabs(x) : dabs(y) )

PLA_DOUBLE_COMPLEX double_complex_inverse( PLA_DOUBLE_COMPLEX x )
{
  PLA_DOUBLE_COMPLEX answer;
  double x2_p_y2, dmax;

  if ( x.real == 0.0 && x.imaginary == 0.0 )
    PLA_Abort( "complex divide by 0", __LINE__, __FILE__ );

  dmax = abs_max( x.real, x.imaginary ); 
 
  x.real = x.real/dmax;
  x.imaginary = x.imaginary/dmax;

  x2_p_y2 = x.real * x.real + x.imaginary * x.imaginary;
  answer.real      =   x.real      / x2_p_y2;
  answer.imaginary = - x.imaginary / x2_p_y2;

  answer.real /= dmax;
  answer.imaginary /= dmax;

  return answer;
}


PLA_COMPLEX complex_inverse( PLA_COMPLEX x )
{
  PLA_COMPLEX answer;
  float x2_p_y2, dmax;

  if ( x.real == 0.0 && x.imaginary == 0.0 )
    PLA_Abort( "complex divide by 0", __LINE__, __FILE__ );

  dmax = abs_max( x.real, x.imaginary );
  x.real = x.real/dmax;
  x.imaginary = x.imaginary/dmax;

  x2_p_y2 = x.real * x.real + x.imaginary * x.imaginary;
  answer.real      =   x.real      / x2_p_y2;
  answer.imaginary = - x.imaginary / x2_p_y2;

  answer.real /= dmax;
  answer.imaginary /= dmax;

  return answer;
}
  
