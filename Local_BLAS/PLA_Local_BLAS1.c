/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

/* Note : I am really unsure what to do if the strides are not */
/*        the same.  Should y's field indicate the new stride? */
/*        Should the function exit?  PLA_BAD_STRIDE?           */

int PLA_Local_copy( PLA_Obj x, PLA_Obj y )
{
   int value = PLA_SUCCESS;
   int original_error = 0;

   int error_value;
   int local_length_x, local_width_x;
   int local_length_y, local_width_y;
   int ld_y_buf, stride_y_buf;
   int rows_in_x_buf, cols_in_x_buf;
   int orientation_x, orientation_y;
   void *y_buf;
   long int size_x, size_y;

   /* Perform parameter checking */

   if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
     value = PLA_Local_copy_enter( x, y );

   PLA_Obj_local_length(x, &local_length_x);
   PLA_Obj_local_width (x, &local_width_x);

   PLA_Obj_local_length(y, &local_length_y);
   PLA_Obj_local_width (y, &local_width_y);

   PLA_Obj_project_onto(x, &orientation_x);
   PLA_Obj_project_onto(y, &orientation_y);
   PLA_Obj_local_buffer(y, &y_buf);
   PLA_Obj_local_stride(y, &stride_y_buf);
   PLA_Obj_local_ldim  (y, &ld_y_buf);

   size_x = (local_length_x * local_width_x);
   size_y = (local_length_y * local_width_y);
   if(  size_x != 0 && size_y != 0 ) {
/*      if(orientation_x == orientation_y)  */
         PLA_Obj_get_local_contents(x, PLA_NO_TRANS,
             &rows_in_x_buf, &cols_in_x_buf, y_buf, ld_y_buf, stride_y_buf);
/*      else
         PLA_Obj_get_local_contents(x, PLA_TRANS,
             &rows_in_x_buf, &cols_in_x_buf, y_buf, ld_y_buf, stride_y_buf); */
    }

/*   PLA_Local_copy_exit(x, y); */

   if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
     value = PLA_Local_copy_exit( x, y );

   return value;
}


int PLA_Local_scal(PLA_Obj alpha, PLA_Obj x)
{
   int value = PLA_SUCCESS;
   int original_error = 0;
   int error_value, i;
   int local_length_x, local_width_x, ld_x;
   int ld_x_buf, stride_x, orientation_x;
   MPI_Datatype datatype_x;
   void *x_buffer, *x_buffer_cur, *alpha_buffer;

   /* Perform parameter checking */

   if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
     value = PLA_Local_scal_enter( alpha, x );

   PLA_Obj_local_length(x, &local_length_x);
   PLA_Obj_local_width (x, &local_width_x);

   PLA_Obj_datatype(x, &datatype_x);

   PLA_Obj_local_ldim(x, &ld_x);
   PLA_Obj_local_stride(x, &stride_x);

   PLA_Obj_project_onto(x, &orientation_x);

   PLA_Obj_local_buffer(x, &x_buffer);
   x_buffer_cur = x_buffer;

   PLA_Obj_local_buffer(alpha, &alpha_buffer);

   if(local_length_x * local_width_x > 0)
   {
      if( PLA_PROJ_ONTO_ROW == orientation_x )
      {
         if(MPI_DOUBLE == datatype_x)
            for(i = 0; i < local_length_x; i++)
            {
               PLA_dscal(&local_width_x, (double*)alpha_buffer, (double*)x_buffer_cur, &ld_x);
               x_buffer_cur = (void*)(((double*)x_buffer_cur)+stride_x);
            }
         else if(MPI_FLOAT == datatype_x)
            for(i = 0; i < local_length_x; i++)
            {
               PLA_sscal(&local_width_x, (float*)alpha_buffer, (float*)x_buffer_cur, &ld_x);
               x_buffer_cur = (void*)(((float*)x_buffer_cur)+stride_x);
            }
         else if(MPI_COMPLEX == datatype_x)
            for(i = 0; i < local_length_x; i++)
            {
               PLA_cscal(&local_width_x, (PLA_COMPLEX*)alpha_buffer,
                         (PLA_COMPLEX*)x_buffer_cur, &ld_x);
               x_buffer_cur = (void*)(((PLA_COMPLEX*)x_buffer_cur)+stride_x);
            }
         else if(MPI_DOUBLE_COMPLEX == datatype_x)
            for(i = 0; i < local_length_x; i++)
            {
               PLA_zscal(&local_width_x, (PLA_DOUBLE_COMPLEX*)alpha_buffer, 
			 (PLA_DOUBLE_COMPLEX*)x_buffer_cur, &ld_x);
               x_buffer_cur = (void*)(((PLA_DOUBLE_COMPLEX*)x_buffer_cur)+stride_x);
            }

      }
      else
      {
         if(MPI_DOUBLE == datatype_x)
            for(i = 0; i < local_width_x; i++)
            {
               PLA_dscal(&local_length_x, (double*)alpha_buffer, 
			 (double*)x_buffer_cur, &stride_x);
               x_buffer_cur = (void*)(((double*)x_buffer_cur)+ld_x);
            }
         else if(MPI_FLOAT == datatype_x)
            for(i = 0; i < local_width_x; i++)
            {
               PLA_sscal(&local_length_x, (float*)alpha_buffer, 
			 (float*)x_buffer_cur, &stride_x);
               x_buffer_cur = (void*)(((float*)x_buffer_cur)+ld_x);
            }
/* I realize that this MAY not work on the Cray -- not really sure */
         else if(MPI_COMPLEX == datatype_x)
            for(i = 0; i < local_width_x; i++)
            {
               PLA_cscal(&local_length_x, (PLA_COMPLEX*)alpha_buffer,
                         (PLA_COMPLEX*)x_buffer_cur, &stride_x);
               x_buffer_cur = (void*)(((PLA_COMPLEX*)x_buffer_cur)+ld_x); /*complex == 2 floats */
            }
         else if(MPI_DOUBLE_COMPLEX == datatype_x)
            for(i = 0; i < local_width_x; i++)
            {
               PLA_zscal(&local_length_x, (PLA_DOUBLE_COMPLEX*)alpha_buffer,
                         (PLA_DOUBLE_COMPLEX*)x_buffer_cur, &stride_x);
               x_buffer_cur = (void*)(((PLA_DOUBLE_COMPLEX*)x_buffer_cur)+ld_x); /* dcomp == 2 doubles */
            }
      }
   }

   if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
     value = PLA_Local_scal_exit( alpha, x );

   return value;
}

extern PLA_DOUBLE_COMPLEX double_complex_inverse( PLA_DOUBLE_COMPLEX );
extern PLA_COMPLEX               complex_inverse( PLA_COMPLEX );

int PLA_Local_inv_scal(PLA_Obj alpha, PLA_Obj x)
{
   int value = PLA_SUCCESS;
   int original_error = 0;
   int error_value, i;
   int local_length_x, local_width_x, ld_x;
   int ld_x_buf, stride_x, orientation_x;
   MPI_Datatype datatype_x;
   void *x_buffer, *x_buffer_cur, *alpha_buffer;

   /* Perform parameter checking */

   if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
     value = PLA_Local_inv_scal_enter( alpha, x );

   PLA_Obj_local_length(x, &local_length_x);
   PLA_Obj_local_width (x, &local_width_x);

   PLA_Obj_datatype(x, &datatype_x);

   PLA_Obj_local_ldim(x, &ld_x);
   PLA_Obj_local_stride(x, &stride_x);

   PLA_Obj_project_onto(x, &orientation_x);

   PLA_Obj_local_buffer(x, &x_buffer);
   x_buffer_cur = x_buffer;

   PLA_Obj_local_buffer(alpha, &alpha_buffer);

   if(local_length_x * local_width_x > 0)
   {
      if( PLA_PROJ_ONTO_ROW == orientation_x )
      {
	if (MPI_DOUBLE == datatype_x){
	  double alpha_inv;

	  alpha_inv = 1.0/(*((double *) alpha_buffer));
	  for(i = 0; i < local_length_x; i++)
            {
	      PLA_dscal(&local_width_x, &alpha_inv, (double*)x_buffer_cur, &ld_x);
	      x_buffer_cur = (void*)(((double*)x_buffer_cur)+stride_x);
            }
	}
	else if(MPI_FLOAT == datatype_x){
	  float alpha_inv;

	  alpha_inv = 1.0/(*((float *) alpha_buffer));
	  for(i = 0; i < local_length_x; i++)
            {
	      PLA_sscal(&local_width_x, &alpha_inv, (float*)x_buffer_cur, &ld_x);
	      x_buffer_cur = (void*)(((float*)x_buffer_cur)+stride_x);
            }
	}
	else if(MPI_COMPLEX == datatype_x){
	  PLA_COMPLEX alpha_inverse;

	  alpha_inverse = complex_inverse( *( PLA_COMPLEX * ) alpha_buffer );

	  for(i = 0; i < local_length_x; i++)
            {
	      PLA_cscal(&local_width_x, &alpha_inverse,
			(PLA_COMPLEX*)x_buffer_cur, &ld_x);
	      x_buffer_cur = (void*)(((PLA_COMPLEX*)x_buffer_cur)+stride_x);
            }
	}
	else if(MPI_DOUBLE_COMPLEX == datatype_x){
	  PLA_DOUBLE_COMPLEX alpha_inverse;

	  alpha_inverse = double_complex_inverse( 
			   *( PLA_DOUBLE_COMPLEX * ) alpha_buffer );

	  for(i = 0; i < local_length_x; i++)
            {
	      PLA_zscal(&local_width_x, &alpha_inverse, 
			(PLA_DOUBLE_COMPLEX*)x_buffer_cur, &ld_x);
	      x_buffer_cur = (void*)(((PLA_DOUBLE_COMPLEX*)x_buffer_cur)+stride_x);
            }
	}
      }
      else
      {
	if (MPI_DOUBLE == datatype_x){
	  double alpha_inv;

	  alpha_inv = 1.0/(*((double *) alpha_buffer));
	  for(i = 0; i < local_width_x; i++)
            {
	      PLA_dscal(&local_length_x, &alpha_inv,
			(double*)x_buffer_cur, &stride_x);
	      x_buffer_cur = (void*)(((double*)x_buffer_cur)+ld_x);
            }
	}
	else if(MPI_FLOAT == datatype_x){
	  float alpha_inv;

	  alpha_inv = 1.0/(*((float *) alpha_buffer));
	  for(i = 0; i < local_width_x; i++)
            {
	      PLA_sscal(&local_length_x, &alpha_inv,
			(float*)x_buffer_cur, &stride_x);
	      x_buffer_cur = (void*)(((float*)x_buffer_cur)+ld_x);
            }
	}
	else if(MPI_COMPLEX == datatype_x){
	  PLA_COMPLEX alpha_inverse;

	  alpha_inverse = complex_inverse( 
			   *( PLA_COMPLEX * ) alpha_buffer );

	  for(i = 0; i < local_width_x; i++)
            {
	      PLA_cscal(&local_length_x, &alpha_inverse,
			(PLA_COMPLEX*)x_buffer_cur, &stride_x);
	      x_buffer_cur = (void*)(((PLA_COMPLEX*)x_buffer_cur)+ld_x); /*complex == 2 floats */
            }
	}
	else if(MPI_DOUBLE_COMPLEX == datatype_x){
	  PLA_DOUBLE_COMPLEX alpha_inverse;

	  alpha_inverse = double_complex_inverse( 
			   *( PLA_DOUBLE_COMPLEX * ) alpha_buffer );
	  for(i = 0; i < local_width_x; i++)
            {
	      PLA_zscal(&local_length_x, &alpha_inverse,
			(PLA_DOUBLE_COMPLEX*)x_buffer_cur, &stride_x);
	      x_buffer_cur =(void*)(((PLA_DOUBLE_COMPLEX*)x_buffer_cur)+ld_x);
	    }
	}
      }
   }

   if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
     value = PLA_Local_inv_scal_exit( alpha, x );

   return value;
}


int PLA_Local_add(PLA_Obj x, PLA_Obj y)
{
  PLA_Obj
    one = NULL;

  PLA_Create_constants_conf_to( x, NULL, NULL, &one );

  PLA_Local_axpy( one, x, y );

  PLA_Obj_free( &one );

  return PLA_SUCCESS;
}


int PLA_Local_axpy( PLA_Obj alpha, PLA_Obj x, PLA_Obj y)
{
  int 
    local_length, local_width,
    ldim_x, ldim_y,
    typesize, i_one = 1, j;

  char
    *buf_x, *buf_y, *buf_alpha;

  MPI_Datatype 
    datatype;
  
  PLA_Obj_local_length(x, &local_length);
  PLA_Obj_local_width (x, &local_width);
  
  if ( local_length * local_width != 0 ){
    PLA_Obj_datatype( x, &datatype );
    MPI_Type_size( datatype, &typesize );
    PLA_Obj_local_buffer(x, ( void ** ) &buf_x);
    PLA_Obj_local_ldim  (x, &ldim_x );
    PLA_Obj_local_buffer(y, ( void ** ) &buf_y);
    PLA_Obj_local_ldim  (y, &ldim_y );
    PLA_Obj_local_buffer(alpha, ( void ** ) &buf_alpha);

    if ( datatype == MPI_DOUBLE ){
      for ( j=0; j<local_width; j++ ){
	PLA_daxpy( &local_length, ( double * ) buf_alpha, 
		                  ( double * ) buf_x, &i_one, 
		                  ( double * ) buf_y, &i_one );
	buf_x += ldim_x * typesize;
	buf_y += ldim_y * typesize;
      }
    }
    else if ( datatype == MPI_FLOAT ){
      for ( j=0; j<local_width; j++ ){
	PLA_saxpy( &local_length, ( float * ) buf_alpha, 
		                  ( float * ) buf_x, &i_one, 
		                  ( float * ) buf_y, &i_one );
	buf_x += ldim_x * typesize;
	buf_y += ldim_y * typesize;
      }
    }
    else if ( datatype == MPI_DOUBLE_COMPLEX ){
      for ( j=0; j<local_width; j++ ){
	PLA_zaxpy( &local_length, ( PLA_DOUBLE_COMPLEX * ) buf_alpha, 
		                  ( PLA_DOUBLE_COMPLEX * ) buf_x, &i_one, 
		                  ( PLA_DOUBLE_COMPLEX * ) buf_y, &i_one );
	buf_x += ldim_x * typesize;
	buf_y += ldim_y * typesize;
      }
    }
    else if ( datatype == MPI_COMPLEX ){
      for ( j=0; j<local_width; j++ ){
	PLA_caxpy( &local_length, ( PLA_COMPLEX * ) buf_alpha, 
		                  ( PLA_COMPLEX * ) buf_x, &i_one, 
		                  ( PLA_COMPLEX * ) buf_y, &i_one );
	buf_x += ldim_x * typesize;
	buf_y += ldim_y * typesize;
      }
    }
  }

  return PLA_SUCCESS;
}


int PLA_Local_nrm2( PLA_Obj x, PLA_Obj alpha )
{
  int 
    local_length, local_width,
    ldim_x, istride;

  char
    *buf_x, *buf_alpha;

  MPI_Datatype 
    datatype;

  PLA_Obj_set_to_zero( alpha );

  PLA_Obj_local_length( x, &local_length);
  PLA_Obj_local_width ( x, &local_width);
  PLA_Obj_local_ldim  ( x, &ldim_x );
  
  if ( local_length * local_width != 0 ){
    if ( local_length == 1 ){
      istride = ldim_x;
      local_length = local_width;
    }
    else
      istride = 1;

    PLA_Obj_datatype( x, &datatype );
    PLA_Obj_local_buffer(x, ( void ** ) &buf_x);
    PLA_Obj_local_buffer(alpha, ( void ** ) &buf_alpha);

    if ( datatype == MPI_DOUBLE ){
      *( ( double * ) buf_alpha ) = 
	PLA_dnrm2( &local_length, ( double *) buf_x, &istride );
    }
    else if ( datatype == MPI_FLOAT ){
      *( ( float * ) buf_alpha ) = 
	PLA_snrm2( &local_length, ( float *) buf_x, &istride );
    }
    else if ( datatype == MPI_DOUBLE_COMPLEX ){
      *( ( double * ) buf_alpha ) = 
	PLA_dznrm2( &local_length, ( PLA_DOUBLE_COMPLEX *) buf_x, &istride );
    }
    else if ( datatype == MPI_COMPLEX ){
       *( ( float * ) buf_alpha ) = 
	      PLA_scnrm2( &local_length, ( PLA_COMPLEX *) buf_x, &istride ); 
    }
  }

  return PLA_SUCCESS;
}

     

int PLA_Local_dot( PLA_Obj x, PLA_Obj y, PLA_Obj alpha )
{
  int 
    local_length_x, local_width_x,
    local_length_y, local_width_y,
    ldim_x, ldim_y, 
    istride_x, istride_y;

  char
    *buf_x, *buf_y, *buf_alpha;

  MPI_Datatype 
    datatype;

  PLA_Obj_set_to_zero( alpha );

  PLA_Obj_local_length( x, &local_length_x);
  PLA_Obj_local_width ( x, &local_width_x);
  PLA_Obj_local_ldim  ( x, &ldim_x );
  PLA_Obj_local_length( y, &local_length_y);
  PLA_Obj_local_width ( y, &local_width_y);
  PLA_Obj_local_ldim  ( y, &ldim_y );
  
  if ( local_length_x * local_width_x != 0 ){
    if ( local_length_x == 1 ){
      istride_x = ldim_x;
      local_length_x = local_width_x;
    }
    else
      istride_x = 1;

    if ( local_length_y == 1 ){
      istride_y = ldim_y;
      local_length_y = local_width_y;
    }
    else
      istride_y = 1;

    PLA_Obj_datatype( x, &datatype );
    PLA_Obj_local_buffer(x, ( void ** ) &buf_x);
    PLA_Obj_local_buffer(y, ( void ** ) &buf_y);
    PLA_Obj_local_buffer(alpha, ( void ** ) &buf_alpha);

    if ( datatype == MPI_DOUBLE ){
      *( ( double * ) buf_alpha ) = 
	PLA_ddot( &local_length_x, ( double *) buf_x, &istride_x,
                                   ( double *) buf_y, &istride_y );
    }
/*    else if ( datatype == MPI_FLOAT ){
      *( ( float * ) buf_alpha ) = 
	PLA_sdot( &local_length_x, ( float *) buf_x, &istride_x,
                                   ( float *) buf_y, &istride_y );
    } */
    else {
      PLA_Abort( "datatype not yet supported", __LINE__, __FILE__ );
    }
  }

  return PLA_SUCCESS;
}

     

int PLA_Local_asum( PLA_Obj x, PLA_Obj alpha )
{
  int 
    local_length, local_width,
    ldim_x, istride;

  char
    *buf_x, *buf_alpha;

  MPI_Datatype 
    datatype;

  PLA_Obj_set_to_zero( alpha );

  PLA_Obj_local_length( x, &local_length);
  PLA_Obj_local_width ( x, &local_width);
  PLA_Obj_local_ldim  ( x, &ldim_x );
  
  if ( local_length * local_width != 0 ){
    if ( local_length == 1 ){
      istride = ldim_x;
      local_length = local_width;
    }
    else
      istride = 1;

    PLA_Obj_datatype( x, &datatype );
    PLA_Obj_local_buffer(x, ( void ** ) &buf_x);
    PLA_Obj_local_buffer(alpha, ( void ** ) &buf_alpha);

    if ( datatype == MPI_DOUBLE ){
      *( ( double * ) buf_alpha ) = 
	PLA_dasum( &local_length, ( double *) buf_x, &istride );
    }
    else if ( datatype == MPI_FLOAT ){
      *( ( float * ) buf_alpha ) = 
	PLA_sasum( &local_length, ( float *) buf_x, &istride );
    }
    else if ( datatype == MPI_DOUBLE_COMPLEX ){
      *( ( double * ) buf_alpha ) = 
	PLA_dzasum( &local_length, ( PLA_DOUBLE_COMPLEX *) buf_x, &istride );
    }
    else if ( datatype == MPI_COMPLEX ){
       *( ( float * ) buf_alpha ) = 
	      PLA_scasum( &local_length, ( PLA_COMPLEX *) buf_x, &istride ); 
    }
  }

  return PLA_SUCCESS;
}


int PLA_Local_swap( PLA_Obj x, PLA_Obj y)
{
  int 
    local_length_x, local_width_x, local_length_y, local_width_y,
    stride_x, stride_y;

  void
    *buf_x, *buf_y;

  MPI_Datatype 
    datatype;
  
  PLA_Obj_local_length(x, &local_length_x);
  PLA_Obj_local_width (x, &local_width_x);
  

  PLA_Obj_local_length(y, &local_length_y);
  PLA_Obj_local_width (y, &local_width_y);

  if ( local_length_x != 0 && local_width_x != 0 &&
       local_length_y != 0 && local_width_y != 0 ){
    PLA_Obj_datatype( x, &datatype );

    PLA_Obj_local_buffer(x, &buf_x);

    if ( local_length_x == 1 ){
      PLA_Obj_local_ldim  (x, &stride_x );
      local_length_x = local_width_x;
    }
    else
      stride_x = 1;

    PLA_Obj_local_buffer(y, &buf_y);

    if ( local_length_y == 1 )
      PLA_Obj_local_ldim  (y, &stride_y );
    else
      stride_y = 1;

    if ( datatype == MPI_DOUBLE ){
      PLA_dswap( &local_length_x, ( double * ) buf_x, &stride_x,
		                   ( double * ) buf_y, &stride_y );
    }
    else if ( datatype == MPI_FLOAT ){
      PLA_sswap( &local_length_x, ( float * ) buf_x, &stride_x,
		                   ( float * ) buf_y, &stride_y );
    }
    else if ( datatype == MPI_DOUBLE_COMPLEX ){
      PLA_zswap( &local_length_x, ( PLA_DOUBLE_COMPLEX * ) buf_x, &stride_x,
     	                   ( PLA_DOUBLE_COMPLEX * ) buf_y, &stride_y );
    }
    else if ( datatype == MPI_COMPLEX ){
      PLA_cswap( &local_length_x, ( PLA_COMPLEX * ) buf_x, &stride_x,
     	                   ( PLA_COMPLEX * ) buf_y, &stride_y );
    }
    else if ( datatype == MPI_INT ){
      int 
	temp, i, *x_p, *y_p;
      
      x_p = (int *) buf_x;
      y_p = (int *) buf_y;

      for ( i=0; i<local_length_x; i++ ){
	temp = *x_p;
	*x_p = *y_p;
	*y_p = temp;
	
	x_p += stride_x;
	y_p += stride_y;
      }
    }
  }

  return PLA_SUCCESS;
}
