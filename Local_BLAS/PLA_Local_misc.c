/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"
#include "math.h"

int PLA_Local_sign( PLA_Obj x, PLA_Obj sign_x )
{
  int local_length, local_width;
  void *buff_x, *buff_sign_x;
  MPI_Datatype datatype;

  PLA_Obj_local_length( x, &local_length );
  PLA_Obj_local_width(  x, &local_width );

  if ( local_length == 1 && local_width == 1 ){
    PLA_Obj_datatype( x, &datatype );
    PLA_Obj_local_buffer( x, &buff_x );
    PLA_Obj_local_buffer( sign_x, &buff_sign_x );

    if ( datatype == MPI_DOUBLE )
      if ( *( (double *) buff_x ) < 0.0 )
	PLA_Obj_set_to_minus_one( sign_x );
      else
	PLA_Obj_set_to_one( sign_x );

    else if ( datatype == MPI_FLOAT)
      if ( *( (float *) buff_x ) < 0.0 )
	PLA_Obj_set_to_minus_one( sign_x );
      else
	PLA_Obj_set_to_one( sign_x );

    else if ( datatype == MPI_COMPLEX ){
      float re, im, nrm;

      re = *(float *)buff_x;
      im = *(((float *)buff_x) + 1);

      nrm = sqrt(re*re+im*im);

      if ( nrm != 0.0 ){
	*(float *)buff_sign_x = re/nrm;
	*(((float *)buff_sign_x) + 1) = im/nrm;
      }
      else
	PLA_Obj_set_to_one ( sign_x);
      
    }

    else if ( datatype == MPI_DOUBLE_COMPLEX ){
      double re, im, nrm;

      re = *(double *)buff_x;
      im = *(((double *)buff_x) + 1);

      nrm = sqrt(re*re+im*im);

      if ( nrm != 0.0 ){
	*(double *)buff_sign_x = re/nrm;
	*(((double *)buff_sign_x) + 1) = im/nrm;
      }
      else
	PLA_Obj_set_to_one ( sign_x);
      
    }

    else {
      printf("PLA_Local_sign: datatype not implemented\n");
      exit (0 );
    }
  }
  else
    PLA_Obj_set_to_zero( sign_x );

  return PLA_SUCCESS;
}
      

int PLA_Local_invert_sign( PLA_Obj alpha )
{
  MPI_Datatype datatype;
  int result=PLA_SUCCESS, local_length, local_width, local_ldim;
  void *local_buffer;

  PLA_Obj_local_length( alpha, &local_length );
  PLA_Obj_local_width( alpha, &local_width );

  if ( local_length==0 && local_width == 0 ) return result;

  PLA_Obj_local_ldim ( alpha, &local_ldim);
  PLA_Obj_datatype( alpha, &datatype );
  PLA_Obj_local_buffer ( alpha, &local_buffer);

  return pla_local_invert_sign (local_length, 
				local_width, 
				local_buffer, 
				local_ldim, 
				datatype );
}
  

    
/* perform the sign inversion of a matrix of specified datatype */
int pla_local_invert_sign ( 
  int m, int n, void *local_buffer, int ld, MPI_Datatype datatype) {

  int i, j;
  float *fbuf = (float *)local_buffer;
  double *dbuf = (double *)local_buffer;
  int *ibuf = (int *) local_buffer;

  for ( j=0; j<n; j++) {
    for ( i=0; i<m; i++) {
      if ( datatype == MPI_FLOAT ){
	*(fbuf+i+j*ld) *= -1.0;
      }
      else if ( datatype == MPI_DOUBLE ){
	*(dbuf+i+j*ld) *= -1.0;
      }
      else if ( datatype == MPI_INT ){
	*(ibuf+i+j*ld) *= -1.0;
      }
      else if ( datatype == MPI_COMPLEX ){
	*(fbuf+2*(i+j*ld)) *= -1.0;
	*(fbuf+2*(i+j*ld)+1) *= -1.0;
      }
      else if ( datatype == MPI_DOUBLE_COMPLEX ){
	*(dbuf+2*(i+j*ld)) *= -1.0;
	*(dbuf+2*(i+j*ld)+1) *= -1.0;
      }
      else {
	printf("Invalid datatype in PLA_Local_invert_sign ()\n");
	return -1;
      }  
    }
  }
  return PLA_SUCCESS;
}






