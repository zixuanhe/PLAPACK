/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

#if MANUFACTURE == CRAY
#define    PLA_NB_OP_PAN_PAN_DEFAULT       96
#define    PLA_NB_OP_PAN_MAT_DEFAULT       96
#define    PLA_NB_OP_MAT_PAN_DEFAULT       96
#define    PLA_NB_OP_PAN_TRIAN_MAT_DEFAULT 96
#define    PLA_NB_OP_SYM_PAN_PAN_DEFAULT   96
#define    PLA_NB_OP_TRIAN_MAT_PAN_DEFAULT 96
 
#elif MANUFACTURE == INTEL
#define    PLA_NB_OP_PAN_PAN_DEFAULT        128
#define    PLA_NB_OP_PAN_MAT_DEFAULT        128
#define    PLA_NB_OP_MAT_PAN_DEFAULT        128
#define    PLA_NB_OP_PAN_TRIAN_MAT_DEFAULT  128 
#define    PLA_NB_OP_SYM_PAN_PAN_DEFAULT    128 
#define    PLA_NB_OP_TRIAN_MAT_PAN_DEFAULT  128 

#elif MANUFACTURE == IBM
#define    PLA_NB_OP_PAN_PAN_DEFAULT       100
#define    PLA_NB_OP_PAN_MAT_DEFAULT       100
#define    PLA_NB_OP_MAT_PAN_DEFAULT       100
#define    PLA_NB_OP_PAN_TRIAN_MAT_DEFAULT 100 
#define    PLA_NB_OP_SYM_PAN_PAN_DEFAULT   100 
#define    PLA_NB_OP_TRIAN_MAT_PAN_DEFAULT 100 

#else
#define    PLA_NB_OP_PAN_PAN_DEFAULT        128
#define    PLA_NB_OP_PAN_MAT_DEFAULT        128
#define    PLA_NB_OP_MAT_PAN_DEFAULT        128
#define    PLA_NB_OP_PAN_TRIAN_MAT_DEFAULT  128 
#define    PLA_NB_OP_SYM_PAN_PAN_DEFAULT    128 
#define    PLA_NB_OP_TRIAN_MAT_PAN_DEFAULT  128       
#endif

static int 
  pla2_nb_op_pan_pan_default       = PLA_NB_OP_PAN_PAN_DEFAULT,
  pla2_nb_op_pan_mat_default       = PLA_NB_OP_PAN_MAT_DEFAULT,
  pla2_nb_op_mat_pan_default       = PLA_NB_OP_MAT_PAN_DEFAULT,
  pla2_nb_op_pan_trian_mat_default = PLA_NB_OP_PAN_TRIAN_MAT_DEFAULT,
  pla2_nb_op_sym_pan_pan_default   = PLA_NB_OP_SYM_PAN_PAN_DEFAULT,
  pla2_nb_op_trian_mat_pan_default = PLA_NB_OP_TRIAN_MAT_PAN_DEFAULT;

int PLA_Environ_nb_alg ( int           type, 
			  PLA_Template  templ,
			  int           *nb_alg )
{
  int value = 0;

  switch ( type ) {
  case PLA_OP_PAN_PAN:
    *nb_alg = pla2_nb_op_pan_pan_default;
    break;
  case PLA_OP_PAN_MAT:
    *nb_alg = pla2_nb_op_pan_mat_default;
    break;
  case PLA_OP_MAT_PAN:
    *nb_alg = pla2_nb_op_mat_pan_default;
    break;
  case PLA_OP_PAN_TRIAN_MAT:
    *nb_alg = pla2_nb_op_pan_trian_mat_default;
    break;
  case PLA_OP_SYM_PAN_PAN:
    *nb_alg = pla2_nb_op_sym_pan_pan_default;
    break;
  case PLA_OP_TRIAN_MAT_PAN:
    *nb_alg = pla2_nb_op_trian_mat_pan_default;
    break;
  default:
    PLA_Abort( "Invalid operation type", __LINE__, __FILE__ );
    value = 1;
  }

  return value;
}


int pla_Environ_set_nb_alg ( int type,  
		              int nb_alg )
{
  int value = 0;

  switch ( type ) {
  case PLA_OP_PAN_PAN:
    pla2_nb_op_pan_pan_default = nb_alg;
    break;
  case PLA_OP_PAN_MAT:
    pla2_nb_op_pan_mat_default = nb_alg;
    break;
  case PLA_OP_MAT_PAN:
    pla2_nb_op_mat_pan_default = nb_alg;
    break;
  case PLA_OP_PAN_TRIAN_MAT:
    pla2_nb_op_pan_trian_mat_default = nb_alg;
    break;
  case PLA_OP_SYM_PAN_PAN:
    pla2_nb_op_sym_pan_pan_default = nb_alg;
    break;
  case PLA_OP_TRIAN_MAT_PAN:
    pla2_nb_op_trian_mat_pan_default = nb_alg;
    break;
  case PLA_OP_ALL_ALG:
    pla2_nb_op_pan_pan_default = nb_alg;
    pla2_nb_op_pan_mat_default = nb_alg;
    pla2_nb_op_mat_pan_default = nb_alg;
    pla2_nb_op_pan_trian_mat_default = nb_alg;
    pla2_nb_op_sym_pan_pan_default = nb_alg;
    pla2_nb_op_trian_mat_pan_default = nb_alg;
    break;
  default:
    PLA_Abort( "Invalid operation type", __LINE__, __FILE__ );
    value = 1;
  }


  return value;
}
