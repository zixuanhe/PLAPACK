/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Conjugate( PLA_Obj x )
{
	int original_error = 0;


	int i;
	int object_type_x;
	int local_length_x, local_width_x, ld_x;
	int stride_x;

	MPI_Datatype datatype_x;
        PLA_Template platemplate;
	PLA_Obj Neg_one = NULL;

	void *x_buffer, *x_buffer_cur;
	void *Neg_one_buffer;

	PLA_Obj_local_length(x, &local_length_x);
	PLA_Obj_local_width (x, &local_width_x);
	PLA_Obj_datatype(x, &datatype_x);
        PLA_Obj_template(x, &platemplate);
	if( (MPI_DOUBLE == datatype_x) || (MPI_DOUBLE_COMPLEX == datatype_x))
		PLA_Mscalar_create( MPI_DOUBLE, PLA_ALL_ROWS, PLA_ALL_COLS,
								  1, 1, platemplate, &Neg_one );
	else PLA_Mscalar_create( MPI_FLOAT, PLA_ALL_ROWS, PLA_ALL_COLS,
								  1, 1, platemplate, &Neg_one );
	PLA_Obj_set_to_minus_one(Neg_one);
	PLA_Obj_local_buffer(Neg_one, &Neg_one_buffer);
	PLA_Obj_local_ldim(x, &ld_x);
	PLA_Obj_objtype(x,  &object_type_x);
	PLA_Obj_local_buffer(x, &x_buffer);
		  PLA_Obj_local_stride(x, &stride_x);
	x_buffer_cur = x_buffer;
	if ( (MPI_DOUBLE == datatype_x) || (MPI_FLOAT == datatype_x))
		local_length_x = local_length_x/2; /* treat AS complex -- not negate */
	if(local_length_x * local_width_x > 0)
	{
		 if( (MPI_DOUBLE_COMPLEX == datatype_x) || (MPI_DOUBLE == datatype_x))
		 {
			 stride_x = stride_x * 2;
			 for(i = 0; i < local_width_x; i++)
			 {
				 PLA_dscal(&local_length_x,(double*)Neg_one_buffer,
							  (((double*)(x_buffer_cur))+1),&stride_x);
				 x_buffer_cur = (void*)(((double*)(x_buffer_cur))+ld_x*2);
			 }
		  }
		 if( (MPI_COMPLEX == datatype_x) || (MPI_FLOAT == datatype_x))
		 {
			 stride_x = stride_x * 2;
			 for(i = 0; i < local_width_x; i++)
			 {
				 PLA_sscal(&local_length_x,(float*)Neg_one_buffer,
							  (((float*)(x_buffer_cur))+1),&stride_x);
				 x_buffer_cur = (void*)(((float*)(x_buffer_cur))+ld_x*2);
			 }
		  }
	}
	PLA_Obj_free(&Neg_one);
	return(PLA_SUCCESS);
}
