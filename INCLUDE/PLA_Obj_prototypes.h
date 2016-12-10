/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

int PLA_Obj_attach_buffer 
                  ( PLA_Obj, void *, int, int );
int PLA_Obj_attach_buffer_enter 
                  ( PLA_Obj, void *, int, int );
int PLA_Obj_attach_buffer_exit 
                  ( PLA_Obj, void *, int, int );

/*----------------------------------------------------------------------------

Purpose : Attach a given data buffer to the object.

IN     obj               object to which buffer is to be attached
IN     buffer            buffer to be used for data
IN     ldim              leading dimension of buffer
IN     user_buffer       indicates whether PLAPACK provided buffer,
                         or user provides buffer.

----------------------------------------------------------------------------*/

int PLA_Vector_create   
                  ( MPI_Datatype, int, PLA_Template, int, PLA_Obj * );
int PLA_Vector_create_enter   
                  ( MPI_Datatype, int, PLA_Template, int, PLA_Obj * );
int PLA_Vector_create_exit   
                  ( MPI_Datatype, int, PLA_Template, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create distributed vector.


IN     datatype          datatype of object
IN     global_length     global length of vector
IN     template          template for vector and matrix distribution
IN     global_align      alignment to template
OUT    new_obj           object describing created vector

----------------------------------------------------------------------------*/

int PLA_Mvector_create_without_buffer 
                 ( MPI_Datatype, int, int, PLA_Template, int, PLA_Obj * );
int PLA_Mvector_create_without_buffer_enter
                 ( MPI_Datatype, int, int, PLA_Template, int, PLA_Obj * );
int PLA_Mvector_create_without_buffer_exit
                 ( MPI_Datatype, int, int, PLA_Template, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create distributed multi-vector, without buffer for data.

IN     datatype          datatype of object
IN     global_length     global length of multi-vector
IN     global_width      global width of multi-vector
IN     templ             template for vector and matrix distribution
IN     global_align      alignment to template
OUT    new_obj           object describing created multi-vector

----------------------------------------------------------------------------*/

int PLA_Mvector_create   
           ( MPI_Datatype, int, int, PLA_Template, int, PLA_Obj * );
int PLA_Mvector_create_enter 
           ( MPI_Datatype, int, int, PLA_Template, int, PLA_Obj * );
int PLA_Mvector_create_exit  
           ( MPI_Datatype, int, int, PLA_Template, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create distributed multi-vector.

IN     datatype          datatype of object
IN     global_length     global length of multi-vector
IN     global_width      global width of multi-vector
IN     templ             template for vector and matrix distribution
IN     global_align      alignment to template
OUT    new_obj           object describing created multi-vector

----------------------------------------------------------------------------*/

int PLA_Matrix_create_without_buffer   
           ( MPI_Datatype, int, int, PLA_Template, int, int, PLA_Obj * );
int PLA_Matrix_create_without_buffer_enter   
           ( MPI_Datatype, int, int, PLA_Template, int, int, PLA_Obj * );
int PLA_Matrix_create_without_buffer_exit   
           ( MPI_Datatype, int, int, PLA_Template, int, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create distributed matrix without buffer.

IN     datatype          datatype of object
IN     global_length     global length of matrix
IN     global_width      global width of matrix
IN     template          template for vector and matrix distribution
IN     global_align_row  row alignment to template
IN     global_align_col  column alignment to template
OUT    new_matrix        object describing created matrix

----------------------------------------------------------------------------*/

int PLA_Matrix_create   
          ( MPI_Datatype, int, int, PLA_Template, int, int, PLA_Obj * );
int PLA_Matrix_create_enter   
          ( MPI_Datatype, int, int, PLA_Template, int, int, PLA_Obj * );
int PLA_Matrix_create_exit   
          ( MPI_Datatype, int, int, PLA_Template, int, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create distributed matrix.

IN     datatype          datatype of object
IN     global_length     global length of matrix
IN     global_width      global width of matrix
IN     template          template for vector and matrix distribution
IN     global_align_row  row alignment to template
IN     global_align_col  column alignment to template
OUT    new_matrix        object describing created matrix

----------------------------------------------------------------------------*/

int PLA_Mscalar_create_without_buffer  
          ( MPI_Datatype, int, int, int, int, PLA_Template, PLA_Obj * );
int PLA_Mscalar_create_without_buffer_enter  
          ( MPI_Datatype, int, int, int, int, PLA_Template, PLA_Obj * );
int PLA_Mscalar_create_without_buffer_exit  
          ( MPI_Datatype, int, int, int, int, PLA_Template, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create distributed multiscalar without buffer.


IN     datatype          datatype of object
IN     owner_row         index of row of nodes with contains owning node(s)
IN     owner_col         index of column of nodes with contains owning node(s)
IN     length            length of mscalar
IN     width             width of mscalar
IN     template          template for vector and matrix distribution
OUT    new_mscalar       object describing created multiscalar

----------------------------------------------------------------------------*/

int PLA_Mscalar_create  ( MPI_Datatype, int, int, int, int, PLA_Template, 
			   PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create distributed mutltiscalar.


IN     datatype          datatype of object
IN     owner_row         index of row of nodes with contains owning node(s)
IN     owner_col         index of column of nodes with contains owning node(s)
IN     length            length of mscalar
IN     width             width of mscalar
IN     template          template for vector and matrix distribution
OUT    new_mscalar       object describing created multiscalar

----------------------------------------------------------------------------*/

int PLA_Pvector_create  
          ( MPI_Datatype, int, int, int, PLA_Template, int, PLA_Obj * );
int PLA_Pvector_create_enter  
          ( MPI_Datatype, int, int, int, PLA_Template, int, PLA_Obj * );
int PLA_Pvector_create_exit  
          ( MPI_Datatype, int, int, int, PLA_Template, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create projected vector.


IN     datatype          datatype of object
IN     project_onto      direction onto which to project
IN     owner             index of row or column of nodes wich containts projected vector
IN     global_length     length of (projected) vector
IN     global_align      alignment to template
IN     template          template for vector and matrix distribution
OUT    new_proj_vector   object describing created projected vector

----------------------------------------------------------------------------*/

int PLA_Pmvector_create_without_buffer  
          ( MPI_Datatype, int, int, int, int,
					   PLA_Template, int, PLA_Obj * );
int PLA_Pmvector_create_without_buffer_enter  
          ( MPI_Datatype, int, int, int, int,
					   PLA_Template, int, PLA_Obj * );
int PLA_Pmvector_create_without_buffer_exit  
          ( MPI_Datatype, int, int, int, int,
					   PLA_Template, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create projected multi-vector without buffer.

IN     datatype              datatype of object
IN     project_onto          direction onto which to project
IN     owner                 index of row or column of nodes wich containts projected vector
IN     global_proj_length    length of (projected) vector
IN     global_proj_width     width of (projected) vector
IN     global_align          alignment to template
IN     template              template for vector and matrix distribution
OUT    new_proj_mvector      object describing created projected multi-vector

----------------------------------------------------------------------------*/

int PLA_Pmvector_create  
          ( MPI_Datatype, int, int, int, int, PLA_Template, int, PLA_Obj * );
int PLA_Pmvector_create_enter  
          ( MPI_Datatype, int, int, int, int, PLA_Template, int, PLA_Obj * );
int PLA_Pmvector_create_exit  
          ( MPI_Datatype, int, int, int, int, PLA_Template, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create projected multi-vector.

IN     datatype              datatype of object
IN     project_onto          direction onto which to project
IN     owner                 index of row or column of nodes wich containts projected vector
IN     global_proj_length    length of (projected) vector
IN     global_proj_width     width of (projected) vector
IN     global_align          alignment to template
IN     template              template for vector and matrix distribution
OUT    new_proj_mvector      object describing created projected multi-vector

----------------------------------------------------------------------------*/

int PLA_Obj_free   
          ( PLA_Obj * );
int PLA_Obj_free_enter   
          ( PLA_Obj * );
int PLA_Obj_free_exit   
          ( PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Linear algebra object destructor.


IN/OUT     object      object to be freed
-----------------------------------------------------------------------------*/

int PLA_Obj_objtype  
          ( PLA_Obj, int * );
int PLA_Obj_objtype_enter  
          ( PLA_Obj, int * );
int PLA_Obj_objtype_exit  
          ( PLA_Obj, int * );

/*----------------------------------------------------------------------------

Purpose : Extract object type from object.


IN         obj        object to be queried
OUT        objtype    object type of object

----------------------------------------------------------------------------*/

int PLA_Obj_datatype  
          ( PLA_Obj, MPI_Datatype * );
int PLA_Obj_datatype_enter  
          ( PLA_Obj, MPI_Datatype * );
int PLA_Obj_datatype_exit  
          ( PLA_Obj, MPI_Datatype * );

/*----------------------------------------------------------------------------

Purpose : Extract datatype from object.


IN         obj         object to be queried
OUT        datatype    datatype of object

----------------------------------------------------------------------------*/

int PLA_Obj_template   
          ( PLA_Obj, PLA_Template * );
int PLA_Obj_template_enter   
          ( PLA_Obj, PLA_Template * );
int PLA_Obj_template_exit   
          ( PLA_Obj, PLA_Template * );

/*----------------------------------------------------------------------------

Purpose : Extract template from object.


IN         obj               object to be queried
OUT        templ             template

----------------------------------------------------------------------------*/

int PLA_Obj_global_info   
          ( PLA_Obj, int *, int *, int *, int *, int *, 
			     int *, int * );
int PLA_Obj_global_info_enter   
          ( PLA_Obj, int *, int *, int *, int *, int *, 
			     int *, int * );
int PLA_Obj_global_info_exit   
          ( PLA_Obj, int *, int *, int *, int *, int *, 
			     int *, int * );

/*----------------------------------------------------------------------------

Purpose :  Extract global information from linear algebra object.

Note that the individual fields can
be queried (separately) through calls to PLA_Obj_global(_field-name).


IN         obj                 object to be queried
OUT        global_length       global row dimension of object
OUT        global_width        global column dimension of object
OUT        project_onto        direction of projection
OUT        owner_row           row index of owning node(s)
OUT        owner_col           column index of owning node(s)
OUT        global_align(_row)  (row) alignment to template
OUT        global_align_col    column alignment to template

----------------------------------------------------------------------------*/

int PLA_Obj_global_length 
          ( PLA_Obj, int * );
int PLA_Obj_global_length_enter 
          ( PLA_Obj, int * );
int PLA_Obj_global_length_exit 
          ( PLA_Obj, int * );

/*----------------------------------------------------------------------------

Purpose :  Extract global length from linear algebra object.

IN         obj                 object to be queried
OUT        global_length       global row dimension of object

----------------------------------------------------------------------------*/

int PLA_Obj_global_width 
          ( PLA_Obj, int * );
int PLA_Obj_global_width_enter 
          ( PLA_Obj, int * );
int PLA_Obj_global_width_exit 
          ( PLA_Obj, int * );

/*----------------------------------------------------------------------------

Purpose :  Extract global width from linear algebra object.

IN         obj                 object to be queried
OUT        global_width        global column dimension of object

----------------------------------------------------------------------------*/

int PLA_Obj_project_onto  
          ( PLA_Obj, int * );
int PLA_Obj_project_onto_enter  
          ( PLA_Obj, int * );
int PLA_Obj_project_onto_exit  
          ( PLA_Obj, int * );

/*----------------------------------------------------------------------------

Purpose :  Extract projection information from linear algebra object.


IN         obj                 object to be queried
OUT        project_onto        direction of projection

----------------------------------------------------------------------------*/

int PLA_Obj_owner_row     
          ( PLA_Obj, int * );
int PLA_Obj_owner_row_enter     
          ( PLA_Obj, int * );
int PLA_Obj_owner_row_exit     
          ( PLA_Obj, int * );

/*----------------------------------------------------------------------------

Purpose :  Extract mesh row owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_row           row index of owning node(s)
OUT        owner_col           column index of owning node(s)

----------------------------------------------------------------------------*/

int PLA_Obj_owner_col     
          ( PLA_Obj, int * );
int PLA_Obj_owner_col_enter     
          ( PLA_Obj, int * );
int PLA_Obj_owner_col_exit     
          ( PLA_Obj, int * );

/*----------------------------------------------------------------------------

Purpose :  Extract mesh column owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_col           column index of owning node(s)

----------------------------------------------------------------------------*/

int PLA_Obj_global_align  
          ( PLA_Obj, int * );
int PLA_Obj_global_align_enter  
          ( PLA_Obj, int * );
int PLA_Obj_global_align_exit  
          ( PLA_Obj, int * );

/*----------------------------------------------------------------------------

WARNING: check this out!!

Purpose :  Extract global alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align        Alignment to template

----------------------------------------------------------------------------*/

int PLA_Obj_global_align_row  
          ( PLA_Obj, int * );
int PLA_Obj_global_align_row_enter  
          ( PLA_Obj, int * );
int PLA_Obj_global_align_row_exit  
          ( PLA_Obj, int * );

/*-------------------------------------------------------------------------

Purpose :  Extract global row alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_row  row alignment to template

-------------------------------------------------------------------------*/

int PLA_Obj_global_align_col 
          ( PLA_Obj, int * );
int PLA_Obj_global_align_col_enter 
          ( PLA_Obj, int * );
int PLA_Obj_global_align_col_exit 
          ( PLA_Obj, int * );

/*-------------------------------------------------------------------------

Purpose :  Extract global column alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_col    column alignment to template

-------------------------------------------------------------------------*/

int PLA_Obj_local_info   
          ( PLA_Obj, int *, int *, void **, int *, int * );
int PLA_Obj_local_info_enter   
          ( PLA_Obj, int *, int *, void **, int *, int * );
int PLA_Obj_local_info_exit   
          ( PLA_Obj, int *, int *, void **, int *, int * );

/*----------------------------------------------------------------------------

Purpose :  Extract local information from linear algebra object.

Note that the individual fields can
be queried (separately) through calls to PLA_Obj_local(_field-name).



IN         obj                 object to be queried
OUT        local_length        local row dimension of object
OUT        local_width         local column dimension of object
OUT        local_buffer        address of local data
OUT        local_stride        stride between entries in a column or vector
OUT        local_ldim          leading dimension of array holding local data

----------------------------------------------------------------------------*/

int PLA_Obj_local_length   
          ( PLA_Obj,       int    *local_length );
int PLA_Obj_local_length_enter   
          ( PLA_Obj,       int    *local_length );
int PLA_Obj_local_length_exit   
          ( PLA_Obj,       int    *local_length );

/*----------------------------------------------------------------------------

Purpose :  Extract local length from linear algebra object.

IN         obj                 object to be queried
OUT        local_length        local row dimension of object

----------------------------------------------------------------------------*/

int PLA_Obj_local_width   
          (PLA_Obj   obj,    int       *local_width );
int PLA_Obj_local_width_enter   
          (PLA_Obj   obj,    int       *local_width );
int PLA_Obj_local_width_exit   
          (PLA_Obj   obj,    int       *local_width );

/*----------------------------------------------------------------------------

Purpose :  Extract local width from linear algebra object.

IN         obj                 object to be queried
OUT        local_width         local column dimension of object

----------------------------------------------------------------------------*/

int PLA_Obj_local_buffer  
          (PLA_Obj   obj,    void   **local_buffer );
int PLA_Obj_local_buffer_enter  
          (PLA_Obj   obj,    void   **local_buffer );
int PLA_Obj_local_buffer_exit  
          (PLA_Obj   obj,    void   **local_buffer );

/*----------------------------------------------------------------------------

Purpose :  Extract local buffer from linear algebra object.

IN         obj                 object to be queried
OUT        local_buffer        address of local data

----------------------------------------------------------------------------*/

int PLA_Obj_local_ldim   
          (PLA_Obj   obj,             int    *local_ldim );
int PLA_Obj_local_ldim_enter   
          (PLA_Obj   obj,             int    *local_ldim );
int PLA_Obj_local_ldim_exit   
          (PLA_Obj   obj,             int    *local_ldim );

/*----------------------------------------------------------------------------

Purpose :  Extract local leading dimension from linear algebra object.

IN         obj                 object to be queried
OUT        local_ldim          leading dimension of array holding local data

----------------------------------------------------------------------------*/

int PLA_Obj_local_stride   
          (PLA_Obj   obj,        int    *local_stride );
int PLA_Obj_local_stride_enter   
          (PLA_Obj   obj,        int    *local_stride );
int PLA_Obj_local_stride_exit   
          (PLA_Obj   obj,        int    *local_stride );

/*----------------------------------------------------------------------------

Purpose :  Extract local stride from linear algebra object.

IN         obj                 object to be queried
OUT        local_stride        stride between entries in a column or vector

----------------------------------------------------------------------------*/


int PLA_Obj_get_local_contents   
          ( PLA_Obj, int, int *, int *, 
				    void *, int, int );
int PLA_Obj_get_local_contents_enter   
          ( PLA_Obj, int, int *, int *, 
				    void *, int, int );
int PLA_Obj_get_local_contents_exit   
          ( PLA_Obj, int, int *, int *, 
				    void *, int, int );

/*----------------------------------------------------------------------------

Purpose : Extract the local data from the given object.

IN     laobject          global object
IN     transpose         indicates whether to transpose data
OUT    rows_in_buf       row dimension of extracted data
OUT    cols_in_buf       column dimension of extracted data
OUT    buf               address where data is to be put
IN     leading_dim       leading dimension of buffer where data is put
IN     stride_buf        stride of buffer where data is put

----------------------------------------------------------------------------*/

int PLA_Obj_set_local_contents   
          (int, int, int, void *, int, int, PLA_Obj );
int PLA_Obj_set_local_contents_enter   
          (int, int, int, void *, int, int, PLA_Obj );
int PLA_Obj_set_local_contents_exit   
          (int, int, int, void *, int, int, PLA_Obj );

/*----------------------------------------------------------------------------

Purpose : Set the local data from the given object.


IN        transpose         indicates whether to transpose data
IN        rows_in_buf       row dimension of data buffer
IN        cols_in_buf       column dimension of data buffer
IN        buf               address of data buffer
IN        leading_dim       leading dimension of data buffer
IN        stride_buf        stride of data buffer where data is put
IN/OUT    laobject          global object

----------------------------------------------------------------------------*/

int PLA_Obj_set   
          ( PLA_Obj, MPI_Datatype, void * );
int PLA_Obj_set_enter   
          ( PLA_Obj, MPI_Datatype, void * );
int PLA_Obj_set_exit   
          ( PLA_Obj, MPI_Datatype, void * );

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to value given
by value. The value is cast to the given MPI_Datatype via the datatype
parameter.


IN/OUT          obj          linear algebra object
IN              datatype     datatype of value
IN             *value        value to be used for initialization

----------------------------------------------------------------------------*/

int PLA_Obj_set_to_zero   
          ( PLA_Obj );
int PLA_Obj_set_to_zero_enter   
          ( PLA_Obj );
int PLA_Obj_set_to_zero_exit   
          ( PLA_Obj );

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to zero (0).
The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE, etc.)
where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/

int PLA_Obj_set_to_one   
          ( PLA_Obj );
int PLA_Obj_set_to_one_enter   
          ( PLA_Obj );
int PLA_Obj_set_to_one_exit   
          ( PLA_Obj );

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to one (1).
The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE, etc.)
where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/

int PLA_Obj_set_to_minus_one   
          ( PLA_Obj );
int PLA_Obj_set_to_minus_one_enter   
          ( PLA_Obj );
int PLA_Obj_set_to_minus_one_exit   
          ( PLA_Obj );

/*----------------------------------------------------------------------------

Purpose : Sets all entries in the linear algebra object to negative
one (-1). The value is cast to the appropriate type (MPI_FLOAT, MPI_DOUBLE,
etc.) where the datatype is the same as that of obj.


IN/OUT          obj          linear algebra object

----------------------------------------------------------------------------*/

int PLA_Obj_view   
          ( PLA_Obj, int, int, int, int, PLA_Obj * );
int PLA_Obj_view_enter   
          ( PLA_Obj, int, int, int, int, PLA_Obj * );
int PLA_Obj_view_exit   
          ( PLA_Obj, int, int, int, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create a view (reference) into an existing linear algebra
object.


IN     old_obj          object into which view is taken
IN     global_length    row dimension of view
IN     global_width     column dimension of view
IN     align_row        row index in old object of upper-left-hand element of view
IN     align_col        column index in old object of upper-left-hand element of view
IN/OUT new_obj          created view

----------------------------------------------------------------------------*/

int PLA_Obj_view_all   
          ( PLA_Obj, PLA_Obj * );
int PLA_Obj_view_all_enter   
          ( PLA_Obj, PLA_Obj * );
int PLA_Obj_view_all_exit   
          ( PLA_Obj, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create a view (reference) into all of existing linear algebra
object.


IN          old_obj          object into which view is taken
IN/OUT      new_obj          created view

----------------------------------------------------------------------------*/

int PLA_Obj_view_swap    
          ( PLA_Obj *, PLA_Obj * );
int PLA_Obj_view_swap_enter    
          ( PLA_Obj *, PLA_Obj * );
int PLA_Obj_view_swap_exit    
          ( PLA_Obj *, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Swaps two views.


IN/OUT      obj1          linear algebra object
IN/OUT      obj2          linear algebra object

----------------------------------------------------------------------------*/

int PLA_Obj_split_4  
          ( PLA_Obj, int, int, PLA_Obj *,  PLA_Obj *, 
		                            PLA_Obj *,  PLA_Obj * );
int PLA_Obj_split_4_enter  
          ( PLA_Obj, int, int, PLA_Obj *,  PLA_Obj *, 
		                            PLA_Obj *,  PLA_Obj * );
int PLA_Obj_split_4_exit  
          ( PLA_Obj, int, int, PLA_Obj *,  PLA_Obj *, 
		                            PLA_Obj *,  PLA_Obj * );
/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into four quadrants

IN         obj               object to be split
IN         length            row dimension of block (as indicated below)
IN         width             column dimension of block (as indicated below)
OUT        upper_left_obj    new object for upper-left block
OUT        upper_right_obj   new object for upper-right block
OUT        lower_left_obj    new object for lower-left block
OUT        lower_right_obj   new object for lower-right block

Here the value of length and width determine the size of one block while
their sign determines which block the length and widht refer to.
If :
length >= 0        width >= 0            upper_left_obj
length <  0        width >= 0            lower_left_obj
length >= 0        width <  0            upper_right_obj
length <  0        width <  0            lower_right_obj

----------------------------------------------------------------------------*/

int PLA_Obj_horz_split_2   
          (  PLA_Obj, int, PLA_Obj *, PLA_Obj * );
int PLA_Obj_horz_split_2_enter   
          (  PLA_Obj, int, PLA_Obj *, PLA_Obj * );
int PLA_Obj_horz_split_2_exit   
          (  PLA_Obj, int, PLA_Obj *, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into top and bottom


IN         obj               object to be split
IN         length            row dimension of block (as indicated below)
OUT        upper_obj         new object for upper block
OUT        lower_obj         new object for lower block

If :
length >= 0     length specifies   upper_obj
length <  0     length specifies   lower_obj

----------------------------------------------------------------------------*/

int PLA_Obj_vert_split_2   
          (  PLA_Obj, int, PLA_Obj *, PLA_Obj * );
int PLA_Obj_vert_split_2_enter   
          (  PLA_Obj, int, PLA_Obj *, PLA_Obj * );
int PLA_Obj_vert_split_2_exit   
          (  PLA_Obj, int, PLA_Obj *, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Split linear algebra object into left and right


IN         obj               object to be split
IN         width             column dimension of block (as indicated below)
OUT        left_obj          new object for left block
OUT        right_obj         new object for right block

If :
width >= 0     width specifies   left_obj
width <  0     width specifies   right_obj

----------------------------------------------------------------------------*/

int PLA_Obj_view_shift   
          ( PLA_Obj, int, int, int, int );
int PLA_Obj_view_shift_enter   
          ( PLA_Obj, int, int, int, int );
int PLA_Obj_view_shift_exit   
          ( PLA_Obj, int, int, int, int );

/*----------------------------------------------------------------------------

Purpose : Shift the boundaries of the view of the linear algebra object



IN/OUT       obj           object to be grown [shifted]
IN           length_top    length to shift top boundary
IN           width_left    width to shift left boundary
IN           width_right   width to shift right boundary
IN           length_bottom length to shift bottom boundary

Note : Positive is down/right negative is up/left (as if you were in quadrant IV
in the cartesian plane.

----------------------------------------------------------------------------*/

int PLA_Obj_split_size   
          ( PLA_Obj, int, int *, int * );
int PLA_Obj_split_size_enter   
          ( PLA_Obj, int, int *, int * );
int PLA_Obj_split_size_exit   
          ( PLA_Obj, int, int *, int * );

/*----------------------------------------------------------------------------

Purpose : Compute size of split to next block boundary.

IN       obj       object to be split
IN       side      side of split
OUT      size      size to next template subblock split
OUT      owner     index of row or column of nodes that owns the split block

----------------------------------------------------------------------------*/

int PLA_Vector_create_conf_to   
          ( PLA_Obj, PLA_Obj *);
int PLA_Vector_create_conf_to_enter   
          ( PLA_Obj, PLA_Obj *);
int PLA_Vector_create_conf_to_exit   
          ( PLA_Obj, PLA_Obj *);

/*----------------------------------------------------------------------------

Purpose : Create vector conformal to given object.


IN      obj          original object
OUT     new_vector   created object

-------------------------------------------------------------------------*/

int PLA_Mvector_create_conf_to   
          ( PLA_Obj, int, PLA_Obj * );
int PLA_Mvector_create_conf_to_enter   
          ( PLA_Obj, int, PLA_Obj * );
int PLA_Mvector_create_conf_to_exit   
          ( PLA_Obj, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create mvector conformal to given object.


IN      obj                original object
IN      global_width       number of vectors in mvector
OUT     new_mvector        created object

-------------------------------------------------------------------------*/

int PLA_Pvector_create_conf_to   
          ( PLA_Obj, int, int, PLA_Obj * );
int PLA_Pvector_create_conf_to_enter   
          ( PLA_Obj, int, int, PLA_Obj * );
int PLA_Pvector_create_conf_to_exit   
          ( PLA_Obj, int, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create projected vector conformal to given object.


IN      obj                original object
IN      project_onto       mesh direction onto which to project
IN      owner              index of row or column of nodes which contains owning nodes
OUT     new_pvector        created object

-------------------------------------------------------------------------*/

int PLA_Pmvector_create_conf_to   
          ( PLA_Obj, int, int, int, PLA_Obj * );
int PLA_Pmvector_create_conf_to_enter   
          ( PLA_Obj, int, int, int, PLA_Obj * );
int PLA_Pmvector_create_conf_to_exit   
          ( PLA_Obj, int, int, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create projected mvector conformal to given object.


IN      obj                original object
IN      project_onto       mesh direction onto which to project
IN      owner              index of row or column of nodes which contains owning nodes
IN      num_vectors        number of vectors in multivector
OUT     new_pmvector       created object

-------------------------------------------------------------------------*/

int PLA_Matrix_create_conf_to   
          ( PLA_Obj, PLA_Obj * );
int PLA_Matrix_create_conf_to_enter   
          ( PLA_Obj, PLA_Obj * );
int PLA_Matrix_create_conf_to_exit   
          ( PLA_Obj, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create matrix conformal to given object (must be a matrix).



IN      matrix       original object
OUT     new_matrix   created object

-------------------------------------------------------------------------*/

int PLA_Mscalar_create_conf_to   
          ( PLA_Obj, int, int, PLA_Obj * );
int PLA_Mscalar_create_conf_to_enter   
          ( PLA_Obj, int, int, PLA_Obj * );
int PLA_Mscalar_create_conf_to_exit   
          ( PLA_Obj, int, int, PLA_Obj * );

/*----------------------------------------------------------------------------

Purpose : Create multiscalar conformal to given object.


IN      obj           original object
IN      owner_row     index of row of nodes that contains owning node(s)
IN      owner_col     index of column of nodes that contains owning node(s)
OUT     new_mscalar   created object (multiscalar)

----------------------------------------------------------------------------*/

int PLA_Obj_objtype_cast   
          ( PLA_Obj, int );
int PLA_Obj_objtype_cast_enter   
          ( PLA_Obj, int );
int PLA_Obj_objtype_cast_exit   
          ( PLA_Obj, int );

/*----------------------------------------------------------------------------

Purpose : Cast the linear algebra object to have given object type.


IN/OUT        obj            original object
IN            objtype        new object type

-------------------------------------------------------------------------*/

int PLA_Obj_set_orientation   
          ( PLA_Obj, int );
int PLA_Obj_set_orientation_enter   
          ( PLA_Obj, int );
int PLA_Obj_set_orientation_exit   
          ( PLA_Obj, int );

/*----------------------------------------------------------------------------

Purpose : Set the orientation of a matrix panel.


IN/OUT      obj              original object
IN          project_onto     annotation indicating orientation

----------------------------------------------------------------------------*/

int PLA_Obj_get_orientation   
          ( PLA_Obj, int * );
int PLA_Obj_get_orientation_enter   
          ( PLA_Obj, int * );
int PLA_Obj_get_orientation_exit   
          ( PLA_Obj, int * );

/*----------------------------------------------------------------------------

Purpose : Return the orientation of a matrix panel.


IN/OUT      obj              original object
IN          project_onto     annotation indicating orientation

----------------------------------------------------------------------------*/

