/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Compute_House_v( PLA_Obj x, PLA_Obj beta )

/*
  Purpose: Compute vector u so that

      ( I + beta u u^T ) x = +- || x ||_2 e_1

  u is normalized to have unit first element and beta = 2 / || u ||_2^2

  Input:  x       --   vector from which to compute Householder transform
                       (PLA_MVECTOR or width 1)
          beta    --   multiscalar in which to store beta 
                       (PLA_MSCALAR of size 1x1 duplicated to all nodes)

  Output: x       --   First element equals +- || x ||_2, the rest is 
                       overwritten with rest of elements of u.
          beta    --   overwritten with beta.

  Return value: PLA_SUCCESS iff transform is successfully computed.
*/
{
  int    x_length, local_length_x1; 
  double d_minus_two = -2.0;
  PLA_Obj x1 = NULL, 
          x_rest = NULL, 
          x1_dup = NULL, 
          sign_x1_dup = NULL,
	  norm_x = NULL, 
          x1_dup_and_norm_x_rest = NULL, 
          norm_x_rest = NULL;
  PLA_Template  templ;
  MPI_Datatype datatype;

  /* Quick return if possible */
  PLA_Obj_global_length(x, &x_length);
  if (x_length == 1) {
    PLA_Obj_set_to_zero( beta );	
    return PLA_SUCCESS;
  }

  PLA_Obj_template( x, &templ );
  PLA_Obj_datatype( x, &datatype );

  /* x1_dup_and_norm_x_rest will equal the vector / x_1         \
                                                  \ norm_x_rest /
     duplicated on all nodes */

  PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS, 
                      2, 1, templ, &x1_dup_and_norm_x_rest );

  /* norm_x will contain the norm of x, updated as x is updated */

  PLA_Mscalar_create( datatype, PLA_ALL_ROWS, PLA_ALL_COLS, 
                      1, 1, templ, &norm_x);

  /* Partition x = / x_1    \
                   \ x_rest /  */
  PLA_Obj_horz_split_2( x, 1,      &x1,
			           &x_rest );

  /* Partition x1_dup_and_norm_x_rest  = / x_1_dup     \
                                         \ norm_x_rest /  */
  PLA_Obj_horz_split_2( x1_dup_and_norm_x_rest, 1, &x1_dup,
		                                   &norm_x_rest );

  /* Copy x1 to a duplicated multiscalar x1_dup */
  PLA_Copy( x1,x1_dup );

  /* Compute norm of x_rest */
  PLA_Nrm2( x_rest, norm_x_rest);

  /* norm of x = norm of / x1          \
                         \ norm_x_rest / */

  PLA_Local_nrm2( x1_dup_and_norm_x_rest, norm_x);

  /* Check if norm_x == 0, and if so, set beta = 0.0 are return */

  if ( PLA_Local_equal_zero ( norm_x ) ){
    PLA_Obj_set_to_zero( beta );
  }
  else{
    /* Let x1_dup = sign( x1_dup ) * norm_x + x1_dup */
    PLA_Mscalar_create_conf_to( x1, PLA_ALL_ROWS, PLA_ALL_COLS, &sign_x1_dup);
    PLA_Local_sign( x1_dup, sign_x1_dup );
    PLA_Local_axpy( sign_x1_dup, norm_x, x1_dup );

  
    /* Let x1 = - sign( x1_dup ) * norm_x */
    PLA_Local_invert_sign( sign_x1_dup );
    PLA_Obj_local_length( x1, &local_length_x1 );
    if ( local_length_x1 == 1 ){
      PLA_Local_copy( norm_x, x1 );
      PLA_Local_scal( sign_x1_dup, x1 );
    }

    /* Let x_rest = x_rest / updated x1_dup, normalizing so that 
       u has a unit first element */
    PLA_Local_inv_scal( x1_dup, x_rest );
    PLA_Local_inv_scal( x1_dup, norm_x_rest );

    /* Compute norm_x = || u ||_2 */
    PLA_Obj_set_to_one( x1_dup );
    PLA_Local_nrm2( x1_dup_and_norm_x_rest, norm_x );  
  
    /* Let beta = - 2 / || u ||_2 / || u ||_2 */
    /*  PLA_Obj_set( beta, MPI_DOUBLE, &d_two ); */
    PLA_Obj_set( beta, MPI_DOUBLE, &d_minus_two ); 
    PLA_Local_inv_scal( norm_x, beta );
    PLA_Local_inv_scal( norm_x, beta );
  }

  /* some clean up */
  PLA_Obj_free( &x1);
  PLA_Obj_free( &x_rest);
  PLA_Obj_free( &x1_dup);
  PLA_Obj_free( &sign_x1_dup);
  PLA_Obj_free( &norm_x);
  PLA_Obj_free( &x1_dup_and_norm_x_rest);
  PLA_Obj_free( &norm_x_rest);

  return PLA_SUCCESS;
}

