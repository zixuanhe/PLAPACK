/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#define PLA_OBJ_OBJTYPE(obj,out) *(out) = (obj)->objtype

/*----------------------------------------------------------------------------

Purpose : Extract object type from object.


IN         obj        object to be queried
OUT        objtype    object type of object

----------------------------------------------------------------------------*/

#define PLA_OBJ_DATATYPE(obj,out) *(out) = (obj)->base_obj->datatype

/*----------------------------------------------------------------------------

Purpose : Extract datatype from object.


IN         obj         object to be queried
OUT        datatype    datatype of object

----------------------------------------------------------------------------*/

#define PLA_OBJ_TEMPLATE(obj,out) *(out) = (obj)->base_obj->templ

/*----------------------------------------------------------------------------

Purpose : Extract template from object.


IN         obj               object to be queried
OUT        templ             template

----------------------------------------------------------------------------*/

/* int PLA_Obj_global_info   ( PLA_Obj, int *, int *, int *, int *, int *, 
			     int *, int * ); */

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

#define PLA_OBJ_GLOBAL_LENGTH(obj,out)*(out) = (obj)->global_length

/*----------------------------------------------------------------------------

Purpose :  Extract global length from linear algebra object.

IN         obj                 object to be queried
OUT        global_length       global row dimension of object

----------------------------------------------------------------------------*/

#define PLA_OBJ_GLOBAL_WIDTH(obj,out) *(out) = (obj)->global_width

/*----------------------------------------------------------------------------

Purpose :  Extract global width from linear algebra object.

IN         obj                 object to be queried
OUT        global_width        global column dimension of object

----------------------------------------------------------------------------*/

#define PLA_OBJ_PROJECT_ONTO(obj,out) *(out) = (obj)->proj_onto

/*----------------------------------------------------------------------------

Purpose :  Extract projection information from linear algebra object.


IN         obj                 object to be queried
OUT        project_onto        direction of projection

----------------------------------------------------------------------------*/

#define PLA_OBJ_OWNER_ROW(obj,out) *(out) = (obj)->owner_row

/*----------------------------------------------------------------------------

Purpose :  Extract mesh row owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_row           row index of owning node(s)
OUT        owner_col           column index of owning node(s)

----------------------------------------------------------------------------*/

#define PLA_OBJ_OWNER_COL(obj,out) *(out) = (obj)->owner_col

/*----------------------------------------------------------------------------

Purpose :  Extract mesh column owner information rom linear algebra object.


IN         obj                 object to be queried
OUT        owner_col           column index of owning node(s)

----------------------------------------------------------------------------*/

#define PLA_OBJ_GLOBAL_ALIGN(obj,out) *(out) = (obj)->global_align_row

/*----------------------------------------------------------------------------

WARNING: check this out!!

Purpose :  Extract global alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align        Alignment to template

----------------------------------------------------------------------------*/

#define PLA_OBJ_GLOBAL_ALIGN_ROW(obj,out) *(out) = (obj)->global_align_row

/*-------------------------------------------------------------------------

Purpose :  Extract global row alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_row  row alignment to template

-------------------------------------------------------------------------*/

#define PLA_OBJ_GLOBAL_ALIGN_COL(obj,out) *(out) = (obj)->global_align_col

/*-------------------------------------------------------------------------

Purpose :  Extract global column alignment from linear algebra object.


IN         obj                 object to be queried
OUT        global_align_col    column alignment to template

-------------------------------------------------------------------------*/

/* int PLA_Obj_local_info ( PLA_Obj, int *, int *, void **, int *, int * );
*/

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

#define PLA_OBJ_LOCAL_LENGTH(obj,out) *(out) = (obj)->local_length

/*----------------------------------------------------------------------------

Purpose :  Extract local length from linear algebra object.

IN         obj                 object to be queried
OUT        local_length        local row dimension of object

----------------------------------------------------------------------------*/

#define PLA_OBJ_LOCAL_WIDTH(obj,out) *(out) = (obj)->local_width

/*----------------------------------------------------------------------------

Purpose :  Extract local width from linear algebra object.

IN         obj                 object to be queried
OUT        local_width         local column dimension of object

----------------------------------------------------------------------------*/

#define PLA_OBJ_LOCAL_BUFFER(obj,out) *(out) = (obj)->local_buffer

/*----------------------------------------------------------------------------

Purpose :  Extract local buffer from linear algebra object.

IN         obj                 object to be queried
OUT        local_buffer        address of local data

----------------------------------------------------------------------------*/

#define PLA_OBJ_LOCAL_LDIM(obj,out) *(out) = (obj)->base_obj->local_ldim

/*----------------------------------------------------------------------------

Purpose :  Extract local leading dimension from linear algebra object.

IN         obj                 object to be queried
OUT        local_ldim          leading dimension of array holding local data

----------------------------------------------------------------------------*/

#define PLA_OBJ_LOCAL_STRIDE(obj,out) *(out) = (obj)->local_stride

/*----------------------------------------------------------------------------

Purpose :  Extract local stride from linear algebra object.

IN         obj                 object to be queried
OUT        local_stride        stride between entries in a column or vector

----------------------------------------------------------------------------*/

#define PLA_OBJ_SIZE_TOP(obj,out) *(out) = (obj)->size_top

/*----------------------------------------------------------------------------

Purpose :  Extract size of top block from linear algebra object.

IN         obj                 object to be queried
OUT        size_top            size of top block

----------------------------------------------------------------------------*/

#define PLA_OBJ_OWNER_TOP(obj,out) *(out) = (obj)->owner_top

/*----------------------------------------------------------------------------

Purpose :  Extract owner of top block from linear algebra object.

IN         obj                 object to be queried
OUT        owner_top           owner of top block

----------------------------------------------------------------------------*/

#define PLA_OBJ_SIZE_BOTTOM(obj,out) *(out) = (obj)->size_bottom

/*----------------------------------------------------------------------------

Purpose :  Extract size of bottom block from linear algebra object.

IN         obj                 object to be queried
OUT        size_bottom         size of bottom block

----------------------------------------------------------------------------*/

#define PLA_OBJ_OWNER_BOTTOM(obj,out) *(out) = (obj)->owner_bottom

/*----------------------------------------------------------------------------

Purpose :  Extract owner of bottom block from linear algebra object.

IN         obj                 object to be queried
OUT        owner_bottom        owner of bottom block

----------------------------------------------------------------------------*/

#define PLA_OBJ_SIZE_LEFT(obj,out) *(out) = (obj)->size_left  

/*----------------------------------------------------------------------------

Purpose :  Extract size of left block from linear algebra object.

IN         obj                 object to be queried
OUT        size_left           size of left block

----------------------------------------------------------------------------*/

#define PLA_OBJ_OWNER_LEFT(obj,out) *(out) = (obj)->owner_left

/*----------------------------------------------------------------------------

Purpose :  Extract owner of left block from linear algebra object.

IN         obj                 object to be queried
OUT        owner_left          owner of left block

----------------------------------------------------------------------------*/

#define PLA_OBJ_SIZE_RIGHT(obj,out) *(out) = (obj)->size_right

/*----------------------------------------------------------------------------

Purpose :  Extract size of right block from linear algebra object.

IN         obj                 object to be queried
OUT        size_right          size of right block

----------------------------------------------------------------------------*/

#define PLA_OBJ_OWNER_RIGHT(obj,out) *(out) = (obj)->owner_right

/*----------------------------------------------------------------------------

Purpose :  Extract owner of right block from linear algebra object.

IN         obj                 object to be queried
OUT        owner_right         owner of right block

----------------------------------------------------------------------------*/

