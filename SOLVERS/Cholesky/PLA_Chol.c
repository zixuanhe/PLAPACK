/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

int PLA_Chol( int uplo, PLA_Obj A )
/*
  Purpose : Parallel Cholesky Factorization

  IN     uplo        integer, PLA_LOWER_TRIANGULAR or PLA_UPPER_TRIANGULAR
  IN/OUT A           matrix to be factored

  Algorith used:

  ******************************************************************
      
      Partition  A = / A_TL ||   *  \
                     | =====  ===== |
                     \ A_BL || A_BR /
     	     where A_TL is 0 x 0 
      while A_BR is not 0 x 0 
         Determine block size b
	 Partition  A = / A_TL ||   *  \    / A_00 ||  *   |  *   \ 
                        | =====  ===== | =  | ====    ====   ==== |
                        \ A_BL || A_BR /    | A_10 || A_11 |  *   |
                                            | ----    ----   ---- |
                                            \ A_20 || A_21 | A_22 /
     	        where A_00 = A_TL and A_11 is b x b 
         Update A_11 <- L_11 = Chol. Fact.( A_11 )
         Update A_21 <- L_21 = A_21 inv( L_11' )
	 Update A_22 <- A_22 - L_21 * L_21'
         Continue with
	            A = / A_TL ||   *  \    / A_00 |  *   ||  *   \ 
                        | =====  ===== | =  | ----   ----    ---- |
                        \ A_BL || A_BR /    | A_10 | A_11 ||  *   |
                                            | ====   ====    ==== |
                                            \ A_20 | A_21 || A_22 /
      endwhile
              
  ******************************************************************

  NOTE:  For details on how to implement parallel Cholesky factorization
         see
	 
	 R. van de Geijn, Using PLAPACK, The MIT Press, 1997.

	 G. Morrow and R. van de Geijn, "Zen and the Art of High-Performance
                 Parallel Computing," http://www.cs.utexas.edu/users/plapack
*/
{
  int 
    value = PLA_SUCCESS,
    nb_alg;

  PLA_Template
    templ;

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Chol_enter( uplo, A );


  /* Get algorithmic block size to be used */
  PLA_Obj_template( A, &templ );
  PLA_Environ_nb_alg( PLA_OP_PAN_PAN, templ, &nb_alg );

  /* Call appropriate subroutine.  The particular one called in tuned to
     perform best for large matrices. */
  value = PLA_Chol_lower_large( nb_alg, A );

  if ( PLA_ERROR_CHECKING )    /* Perform parameter and error checking */
    value = PLA_Chol_exit( uplo, A );

  return value;
}

int PLA_Chol_lower_large( int nb_alg, PLA_Obj A )
/*
  Purpose : Parallel Cholesky Factorization.  This particular version 
            assumes only the lower triangular portion of A is stored and
            is tuned to perform best for large matrices.

  IN     nb_alg      integer, algorithmic block size to be used
  IN/OUT A           matrix to be factored

  Algorith used:

  ******************************************************************
      
      Partition  A = / A_TL ||   *  \
                     | =====  ===== |
                     \ A_BL || A_BR /
     	     where A_TL is 0 x 0 
      while A_BR is not 0 x 0 
         Determine block size b
	 Partition  A = / A_TL ||   *  \    / A_00 ||  *   |  *   \ 
                        | =====  ===== | =  | ====    ====   ==== |
                        \ A_BL || A_BR /    | A_10 || A_11 |  *   |
                                            | ----    ----   ---- |
                                            \ A_20 || A_21 | A_22 /
     	        where A_00 = A_TL and A_11 is b x b 
         Update A_11 <- L_11 = Chol. Fact.( A_11 )
         Update A_21 <- L_21 = A_21 inv( L_11' )
	 Update A_22 <- A_22 - L_21 * L_21'
         Continue with
	            A = / A_TL ||   *  \    / A_00 |  *   ||  *   \ 
                        | =====  ===== | =  | ----   ----    ---- |
                        \ A_BL || A_BR /    | A_10 | A_11 ||  *   |
                                            | ====   ====    ==== |
                                            \ A_20 | A_21 || A_22 /
      endwhile
              
  ******************************************************************

  NOTE:  For details on how to implement parallel Cholesky factorization
         see
	 
	 R. van de Geijn, Using PLAPACK, The MIT Press, 1997.

	 G. Morrow and R. van de Geijn, "Zen and the Art of High-Performance
                 Parallel Computing," http://www.cs.utexas.edu/users/plapack
*/
{
  int       
    value = 0,
    size, k,
    trans;   

  PLA_Obj  
    ABR             = NULL,     
    A11             = NULL,     A21             = NULL,     
    A11_dup         = NULL,     A21_mv         = NULL,
    A21_dpmv        = NULL,     A12_dpmv        = NULL,
    one             = NULL,     minus_one       = NULL;

  MPI_Datatype
    datatype;
  
  PLA_Obj_datatype( A, &datatype );

  if ( ( datatype == MPI_DOUBLE ) || ( datatype == MPI_FLOAT ) )
    trans = PLA_TRANSPOSE;
  else
    trans = PLA_CONJUGATE_TRANSPOSE;

  /* Create usual duplicated scalar constants */
  PLA_Create_constants_conf_to( A, &minus_one, NULL, &one );

  /* View A_BR = A */
  PLA_Obj_view_all( A, &ABR );

  k = 0;
  while ( TRUE ) {
    /* Determine size of current panel */
    PLA_Obj_global_length( ABR, &size );
    if ( 0 == ( size = min( size, nb_alg ) ) ) break;

    /* Partition A_BR = / A_11   *   \
                        \ A_21  A_22 / where A_11 is b x b.
       Notice that A_22 becomes A_BR in the next iteration, so we
       already view this as ABR. */
    PLA_Obj_split_4( ABR, size, size,   &A11, PLA_DUMMY,
                                         &A21, &ABR );

    /* Create a multiscalar so that A_11 can be duplicated to all processors */
    PLA_Mscalar_create_conf_to( A11, PLA_ALL_ROWS, PLA_ALL_COLS, &A11_dup );
    
    /* Duplicate A_11 to all processors */
    PLA_Copy( A11, A11_dup );

    /* Locally on each processor compute the Cholesky factorization of A11 */
    value = PLA_Local_chol( PLA_LOWER_TRIANGULAR, A11_dup );

    /* Check is error was detected in local Cholesky factorization */
    if ( value != 0 ) {
      printf("value = %d\n", value );
      value = k + value;
      break;
    }

    /* Redistribute A_21 so that all processos can participate in updating 
       A_21 <- L_21 = A_21 A_11^T */
    PLA_Mvector_create_conf_to( A21, size, &A21_mv );
    PLA_Copy( A21, A21_mv );

    /* Perform local part of A_21 <- L_21 = A_21 A_11^T */
    PLA_Local_trsm( PLA_SIDE_RIGHT, PLA_LOWER_TRIANGULAR,
		     trans,  PLA_NONUNIT_DIAG,
		     one, A11_dup, A21_mv );

    /* Duplicate A_21 so that the update of A_22 becomes a local operation */
    PLA_Pmvector_create_conf_to( ABR, PLA_PROJ_ONTO_COL, PLA_ALL_COLS, size,
				 &A21_dpmv );
    PLA_Pmvector_create_conf_to( ABR, PLA_PROJ_ONTO_ROW, PLA_ALL_ROWS, size,
				 &A12_dpmv );

    PLA_Copy( A21_mv, A21_dpmv );
    PLA_Copy( A21_mv, A12_dpmv );

    /* Perform local part of the update of A_22 */
    if ( datatype == MPI_DOUBLE || datatype == MPI_FLOAT )
      PLA_Syrk_perform_local_part( PLA_LOWER_TRIANGULAR, 
				    minus_one, A21_dpmv, A12_dpmv,
                                    one,       ABR );
    else{
      PLA_Conjugate( A12_dpmv );
      PLA_Herk_perform_local_part( PLA_LOWER_TRIANGULAR, 
				    minus_one, A21_dpmv, A12_dpmv,
                                    one,       ABR );
    }

    /* Copy A_11 back into place.  Notice: no communication necessary */
    PLA_Copy( A11_dup, A11 );

    /* Copy A_21 back into place.  Notice: no communication necessary */
    PLA_Copy( A21_dpmv, A21 );
    k = k+size;
  }

  /* free temporary objects and views */
  PLA_Obj_free( &ABR );        
  PLA_Obj_free( &A11 );        PLA_Obj_free( &A21 ); 
  PLA_Obj_free( &A11_dup );    PLA_Obj_free( &A21_mv );
  PLA_Obj_free( &A21_dpmv );   PLA_Obj_free( &A12_dpmv );
  PLA_Obj_free( &one );        PLA_Obj_free( &minus_one );
  
  return value;
}

