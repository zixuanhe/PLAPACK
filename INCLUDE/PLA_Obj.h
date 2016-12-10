/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#if MANUFACTURE == CRAY
#include <ffio.h>
#endif

struct PLA_Obj_base_struct{
  MPI_Datatype 
       datatype;           /* Datatype of linear algebra object */
  void 
       *local_buffer;      /* Local buffer for data */
  int
       global_length,      /* Global length of base object */
       global_width,       /* Global width of base object */
       global_align_row,   /* Global alignment to template (row dimension) */ 
       global_align_col,   /* Global alignment to template (col dimension) */ 
       local_ldim,         /* Local leading dimension of buffer */
       user_buffer,        /* Indicates whether buffer was provided by user */
       number_of_views,    /* Number of views into base objects */
       local_buffer_size,  /* Size of buffer attached to object */
       fd,                 /* File descriptor for OOC object */         
       ooc;                /* Indicates if out of core object */
#if MANUFACTURE == CRAY
  struct ffsw 
       ooc_status;         /* Status of async IO with OOC object */
#else
  int  
       ooc_status;
#endif
  
  PLA_Template 
       templ;              /* Template for data distribution */
};
typedef struct PLA_Obj_base_struct *PLA_Obj_base;

struct PLA_Obj_view_struct{
  int 
       objtype,            /* Type of linear algebra object */
       fixed,              /* Indicates if object is used for constants */
       global_length,      /* Global length of linear algebra object */
       global_width,       /* Global width of linear algebra object */
       global_align_row,   /* Global row alignment of TL element */
       global_align_col,   /* Global col alignment of TL element */
       owner_row,          /* Row index of owning node(s) */
       owner_col,          /* Column index of owning node(s) */
       proj_onto,          /* Orientation */

       size_top,           /* Size to next block boundary at top */
       owner_top,          /* Owner of top block */
       size_bottom,        /* Size to next block boundary at bottom */
       owner_bottom,       /* Owner of bottom block */
       size_left,          /* Size to next block boundary at left */
       owner_left,         /* Owner of left block */
       size_right,         /* Size to next block boundary at right */
       owner_right,        /* Owner of right block */

       local_length,       /* Length of local part of object */
       local_width,        /* Width of local part of object */
       local_stride;       /* Stride of data in object */

  void 
       *local_buffer;      /* Local buffer for data */
  
#if DEBUG==1
  PLA_Obj pla_obj;         /* Copy of object using R1.2 */
#endif
  PLA_Obj_base 
       base_obj;           /* Pointer to base object */

  int 
    mode,
    index_in_open_objects;
};
typedef struct PLA_Obj_view_struct *PLA_Obj;

extern int PLA_ERROR_CHECKING;
extern int PLA_CHECK_PARAMETERS;
extern int PLA_CHECK_AGAINST_SEQUENTIAL;
extern int PLA_CHECK_AGAINST_R12;

#define PLA_COMM_ALL         9999
#define PLA_COMM_ROW         9998
#define PLA_COMM_COL         9997
