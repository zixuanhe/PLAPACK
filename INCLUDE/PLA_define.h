/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#define PLA_OBJTYPE_OFFSET 10

#define PLA_VECTOR                 ( PLA_OBJTYPE_OFFSET + 1 )
#define PLA_MVECTOR                ( PLA_OBJTYPE_OFFSET + 2 ) 
#define PLA_MATRIX                 ( PLA_OBJTYPE_OFFSET + 3 ) 
#define PLA_MSCALAR                ( PLA_OBJTYPE_OFFSET + 4 ) 
#define PLA_PVECTOR                ( PLA_OBJTYPE_OFFSET + 5 ) 
#define PLA_PMVECTOR               ( PLA_OBJTYPE_OFFSET + 6 ) 

#define PLA_TRANSPOSE_OFFSET 20

#define PLA_NO_TRANSPOSE           ( PLA_TRANSPOSE_OFFSET + 1 )
#define PLA_NO_TRANS               ( PLA_TRANSPOSE_OFFSET + 1 )
#define PLA_TRANSPOSE              ( PLA_TRANSPOSE_OFFSET + 2 )
#define PLA_TRANS                  ( PLA_TRANSPOSE_OFFSET + 2 )
#define PLA_CONJUGATE_TRANSPOSE    ( PLA_TRANSPOSE_OFFSET + 3 )
#define PLA_CONJ_TRANS             ( PLA_TRANSPOSE_OFFSET + 3 )
/* not in book */
#define PLA_CONJUGATE              ( PLA_TRANSPOSE_OFFSET + 4 ) 
#define PLA_CONJ                   ( PLA_TRANSPOSE_OFFSET + 4 )

#define PLA_OWNER_OFFSET -30

#define PLA_ALL_ROWS               ( PLA_OWNER_OFFSET - 1 )
#define PLA_ALL_COLS               ( PLA_OWNER_OFFSET - 2 )

#define PLA_SIDE_OFFSET 40

#define PLA_SIDE_TOP               ( PLA_SIDE_OFFSET  + 1 )
#define PLA_SIDE_BOTTOM            ( PLA_SIDE_OFFSET  + 2 )
#define PLA_SIDE_LEFT              ( PLA_SIDE_OFFSET  + 3 )
#define PLA_SIDE_RIGHT             ( PLA_SIDE_OFFSET  + 4 )

#define PLA_SHAPE_OFFSET 50

#define PLA_UPPER_TRIANGULAR       ( PLA_SHAPE_OFFSET + 1 )
#define PLA_UP_TRIAN               ( PLA_SHAPE_OFFSET + 1 )
#define PLA_LOWER_TRIANGULAR       ( PLA_SHAPE_OFFSET + 2 )
#define PLA_LOWER_TRIAN            ( PLA_SHAPE_OFFSET + 2 )
#define PLA_SHAPE_GENERAL          ( PLA_SHAPE_OFFSET + 3 )
#define PLA_SHAPE_LOWER_TRAPEZOIDAL ( PLA_SHAPE_OFFSET + 4 )
#define PLA_SHAPE_LOW_TRAP         ( PLA_SHAPE_OFFSET + 4 )
#define PLA_SHAPE_UPPER_TRAPEZOIDAL ( PLA_SHAPE_OFFSET + 5 )
#define PLA_SHAPE_UP_TRAP          ( PLA_SHAPE_OFFSET + 5 )

#define PLA_FLOW_OFFSET 60

#define PLA_DIR_NOFLOW             ( PLA_FLOW_OFFSET + 1 )
#define PLA_DIR_UP                 ( PLA_FLOW_OFFSET + 2 )
#define PLA_DIR_DOWN               ( PLA_FLOW_OFFSET + 3 )
#define PLA_DIR_LEFT               ( PLA_FLOW_OFFSET + 4 )
#define PLA_DIR_RIGHT              ( PLA_FLOW_OFFSET + 5 )
#define PLA_TEMP_COL               ( PLA_FLOW_OFFSET + 6 )
#define PLA_TEMP_ROW               ( PLA_FLOW_OFFSET + 7 )

#define PLA_API_OFFSET 70

#define PLA_MODE_CLOSED            ( PLA_API_OFFSET + 1 )
#define PLA_MODE_OPEN              ( PLA_API_OFFSET + 2 )
#define PLA_STATE_API_ACTIVE       ( PLA_API_OFFSET + 3 )
#define PLA_STATE_NORMAL           ( PLA_API_OFFSET + 4 )

#define PLA_OP_OFFSET 80

#define PLA_OP_MAT_PAN             ( PLA_OP_OFFSET + 1 )
#define PLA_OP_PAN_MAT             ( PLA_OP_OFFSET + 2 )
#define PLA_OP_PAN_PAN             ( PLA_OP_OFFSET + 3 )
#define PLA_OP_PAN_TRIAN_MAT       ( PLA_OP_OFFSET + 4 )
#define PLA_OP_SYM_PAN_PAN         ( PLA_OP_OFFSET + 5 )
#define PLA_OP_TRIAN_MAT_PAN       ( PLA_OP_OFFSET + 6 )
#define PLA_OP_ALL_ALG             ( PLA_OP_OFFSET + 7 )

#define PLA_PROJ_OFFSET 90

#define PLA_PROJ_ONTO_COL          ( PLA_PROJ_OFFSET + 1 )
#define PLA_PROJ_ONTO_ROW          ( PLA_PROJ_OFFSET + 2 )

#define PLA_DIAG_OFFSET 100
#define PLA_UNIT_DIAG              ( PLA_DIAG_OFFSET + 1 )
#define PLA_NONUNIT_DIAG           ( PLA_DIAG_OFFSET + 2 )
#define PLA_ZERO_DIAG              ( PLA_DIAG_OFFSET + 3 )

/* Misc. */

#define PLA_DUMMY         0
#define PLA_SUCCESS       0
#define PLA_ALIGN_FIRST  -200 
#define PLA_DIM_ALL      -201
#define PLA_INHERIT      -202
#define PLA_UNDEFINED    -203


