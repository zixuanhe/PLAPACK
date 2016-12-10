/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#if DEBUG==1
  #include "PLA.h"
#endif

#include "mpi.h"
#include "PLA_Machines.h"
#include "PLA_OOC.h"
#include "PLA_Misc.h"
#include "PLA_Template.h"
#include "PLA_Obj.h"
#include "PLA_define.h"
#include "PLA_Obj_macros.h"
#include "PLA_API_prototypes.h"
#include "PLA_API_util_prototypes.h"
#include "PLA_Diagnostics_prototypes.h"
#include "PLA_Obj_management_prototypes.h"
#include "PLA_Memory_prototypes.h"
#include "PLA_Obj_prototypes.h"
#include "PLA_BLAS.h"
#include "PLA_BLAS_prototypes.h"
#include "PLA_Versions.h"
#include "PLA_Timings.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef NULL
#define NULL 0
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
