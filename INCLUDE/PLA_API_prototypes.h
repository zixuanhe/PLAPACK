/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

/******************************************************************************/

int PLA_API_begin();

/*--------------------------------------------------------------------------

Purpose : Activates the API active state.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_end();

/*--------------------------------------------------------------------------

Purpose : De-activates the API active state.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_state( int    * );

/*--------------------------------------------------------------------------

Purpose : Determines the state of a linear algebra object.

OUT     state    State of linear algebra object

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_Obj_API_open( PLA_Obj ); 

/*--------------------------------------------------------------------------

Purpose : Puts object into asynchronous mode
          (makes PLA_Obj_API_ operations on the object allowable).


IN/OUT     obj    Linear algebra object (to be) placed in asynchronous mode.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_Obj_API_close( PLA_Obj );

/*--------------------------------------------------------------------------

Purpose : Puts object into synchronous mode (disallows PLA_Obj_API_
operations on the object).


IN/OUT     obj    Linear algebra object (to be) taken out of asynchronous mode.

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_Obj_API_mode (PLA_Obj,     int       * );

/*--------------------------------------------------------------------------

Purpose : Returns the mode of the linear algebra object (PLA_MODE_CLOSED
or PLA_MODE_OPEN).


IN     obj       Linear algebra object
OUT    mode      Mode of linear algebra object

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_axpy_vector_to_global(   int, void *, void *, int, PLA_Obj, int );

/*--------------------------------------------------------------------------

Purpose : Add local vector buffer data to global vector object.


IN            size              Length of vector
IN          * alpha             Scalar in axpy operation
IN          * local_obj         Local vector values
IN            local_inc         Local vector increment (stride)
IN/OUT        laobj             Global vector (target linear algebra object)
IN            disp              Global vector displacement

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_axpy_global_to_vector(int, void *, PLA_Obj, int, void *, int );

/*--------------------------------------------------------------------------

Purpose : Add piece of global vector object to local linear algebra
object.


IN            size              Length of vector
IN          * alpha             Scalar in axpy operation
IN            laobj             Global vector object
IN            displ             Displacement  (in global object)
IN/OUT      * local_buf         Local linear algebra object (target) buffer
IN            local_inc         Local vector increment (stride)

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_multi_axpy_vector_to_global( 
		int, int *, void *, void *, int, PLA_Obj, int * );

/*--------------------------------------------------------------------------

Purpose : Add local sub-vector buffer data chunks to global vector
object.


IN            nsub              Number of subvectors
IN          * size              Length of subvectors
IN          * alpha             Scalar in axpy operation
IN          * local_vector      Local vector values
IN            local_stride      Local vector increment (stride)
IN/OUT        laobj             Global vector (target linear algebra object)
IN          * displs            Global vector displacements

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_multi_axpy_global_to_vector(
			int, int *, void *, PLA_Obj, int *, void *, int );

/*--------------------------------------------------------------------------

Purpose : Add pieces of global vector object to local vector.

 IN            nsub              Number of subvectors
IN          * size              Length of subvectors
IN          * alpha             Scalar in axpy operation
IN            laobj             Global vector object
IN          * displs            Displacements  (in global object)
IN/OUT      * local_vector      Local linear algebra object (target) buffer
IN            local_stride      Local vector increment (stride)

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_axpy_matrix_to_global(
		      int, int, void *, void *, int, PLA_Obj, int, int );

/*--------------------------------------------------------------------------

Purpose : Add local matrix buffer data to global matrix object.


IN            m                 Length of matrix
IN            n                 Width of matrix
IN          * alpha             Scalar in axpy operation
IN          * local_buf         Local matrix values (buffer)
IN            local_ldim        Local matrix leading dimension
IN/OUT        A                 Global matrix (target linear algebra object)
IN            disp_row          Global matrix row displacement
IN            disp_col          Global matrix column displacement

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_axpy_global_to_matrix(
		  int, int, void *, PLA_Obj, int, int, void *, int );

/*--------------------------------------------------------------------------

Purpose : Add piece of global matrix object to local linear algebra
object (matrix).


IN            m                 Length of matrix
IN            n                 Width of matrix
IN          * alpha             Scalar in axpy operation
IN            A                 Global object
IN            displ_row         Row displacement  (in global object)
IN            displ_col         Column displacement  (in global object)
IN/OUT      * local_buf      Local linear algebra object (target) buffer
IN            local_ldim         Local matrix leading dimension

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_API_multi_axpy_matrix_to_global(   
	int, int, int *, int *, void *, void *, int, PLA_Obj, int *, int * );

/*--------------------------------------------------------------------------

Purpose : Add local sub-matrix buffer data chunks to global matrix
object.


IN            nsub_row          Number of row subblocks to map
IN            nsub_row          Number of column subblocks to map
IN          * size_row          Row block lengths
IN          * size_col          Column block widths
IN          * alpha             Scalar in axpy operation
IN          * local_matrix      Local matrix values
IN            local_ldim        Local matrix leading dimension
IN/OUT        obj               Global matrix (target linear algebra object)
IN          * disp_row          Global matrix row displacements
IN          * disp_col          Global matrix column displacements

----------------------------------------------------------------------------*/

/****************************************************************************/

int PLA_API_multi_axpy_global_to_matrix(   
	int, int, int *, int *, void *, PLA_Obj, int *, int *, void *, int );

/*--------------------------------------------------------------------------

Purpose : Add pieces of global matrix object to local matrix.



IN            nsub_row          Number of row subblocks to map
IN            nsub_row          Number of column subblocks to map
IN          * size_row          Row block lengths
IN          * size_col          Column block widths
IN          * alpha             Scalar in axpy operation
IN            obj               Global matrix
IN          * disp_row          Global matrix row displacements
IN          * disp_col          Global matrix column displacements
IN/OUT      * local_matrix      Local matrix values
IN            local_ldim        Local matrix leading dimension

----------------------------------------------------------------------------*/

/******************************************************************************/

int PLA_Obj_API_sync( PLA_Obj );

/*--------------------------------------------------------------------------

Purpose : Synchronizes all nodes w.r.t. the given object and completes
all pending "API_axpy" operations on the object.


IN/OUT    laobj      Linear algebra object

----------------------------------------------------------------------------*/

