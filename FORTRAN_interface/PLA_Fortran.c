/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

#include "PLA.h"

#if MANUFACTURE == IBM
#define MPI_Fint int
#define MPI_Comm_c2f(obj) (obj)
#define MPI_Comm_f2c(obj) (obj)
#define MPI_Op_f2c(op)    (op)
#endif


#if MANUFACTURE == CRAY
#define MPI_Fint int
#define MPI_Comm_c2f(obj) (obj)
#define MPI_Comm_f2c(obj) (obj)
#define MPI_Op_f2c(op)    (op)
#endif

#if MANUFACTURE != SGI
#define FINT int
#endif

#if MANUFACTURE == PC
#undef PLA_FOR2C
#define PLA_FOR2C( name ) name ## __ 
#endif

/*************************************************************
   Chapter 2
**************************************************************/

/* 
   Section 2.1: Initializing PLAPACK
*/
 
#if MANUFACTURE == CRAY
void PLA_INIT_F
#else
void PLA_FOR2C( pla_init_f )
#endif
( MPI_Fint *comm, FINT *ierror )
{
  MPI_Comm 
    comm_C = MPI_COMM_NULL; 

  comm_C = MPI_Comm_f2c( *comm );
  *ierror = ( FINT) PLA_Init( comm_C );

  return;
}

#if MANUFACTURE == CRAY
void PLA_FINALIZE_F
#else
void PLA_FOR2C( pla_finalize_f )
#endif
( FINT *ierror )
{
  *ierror = ( FINT ) PLA_Finalize();

  return;
}

#if MANUFACTURE == CRAY
void PLA_INITIALIZED_F
#else
void PLA_FOR2C( pla_initialized_f )
#endif
( FINT *initialized, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Initialized( ( int ) initialized );

  return;
}

#if MANUFACTURE == CRAY
void PLA_BASECOMM_F
#else
void PLA_FOR2C( pla_basecomm_f )
#endif
( MPI_Fint *comm, FINT *ierror )
{
  MPI_Comm 
    comm_C = MPI_COMM_NULL;
  int 
    value;

  value = PLA_Base_comm( &comm_C );
  
  *comm = MPI_Comm_c2f( comm_C );
  
  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_COMM1DTO2D_F
#else
void PLA_FOR2C( pla_comm1dto2d_f )
#endif
( MPI_Fint *comm_in,  FINT *nprows, FINT *npcols, 
			   MPI_Fint *comm_out, FINT *ierror )
{
  MPI_Comm 
    comm_in_C = MPI_COMM_NULL, 
    comm_out_C = MPI_COMM_NULL;
  int value;

  comm_in_C = MPI_Comm_f2c( *comm_in );
  value = PLA_Comm_1D_to_2D( comm_in_C, ( int ) *nprows, ( int ) *npcols, 
			     &comm_out_C );

  *comm_out = MPI_Comm_c2f( comm_out_C );

  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_COMM1DTO2DRATIO_F
#else
void PLA_FOR2C( pla_comm1dto2dratio_f )
#endif
       ( MPI_Fint *comm_in, double *ratio, MPI_Fint *comm_out, FINT *ierror )
{
  MPI_Comm 
    comm_in_C = MPI_COMM_NULL, 
    comm_out_C = MPI_COMM_NULL;
  int 
    value;

  comm_in_C = MPI_Comm_f2c( *comm_in );
  value = PLA_Comm_1D_to_2D_ratio
       ( comm_in_C, ( float ) *ratio, &comm_out_C );

  *comm_out = MPI_Comm_c2f( comm_out_C );

  *ierror = ( FINT ) value;

  return;
}

/* 
   Section 2.2: Distribution Templates
*/

#if MANUFACTURE == CRAY
void PLA_TEMPCREATE_F
#else
void PLA_FOR2C( pla_tempcreate_f )
#endif
( FINT *nb_distr, FINT *zero_or_one, 
			FINT *templ,    FINT *ierror )
{
  *ierror = ( FINT ) PLA_Temp_create( ( int ) *nb_distr, 
				      ( int ) *zero_or_one, 
				      ( PLA_Template *) templ );

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPFREE_F
#else
void PLA_FOR2C( pla_tempfree_f )
#endif
( FINT *templ, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Temp_free( (PLA_Template *) templ );

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMALLINFO_F
#else
void PLA_FOR2C( pla_tempcommallinfo_f )
#endif
( FINT *templ, 
			       MPI_Fint *comm, 
			       FINT *rank, 
			       FINT *numnodes, 
			       FINT *ierror )
{
  MPI_Comm 
    comm_C = MPI_COMM_NULL;
  int 
    value,
    rank_C, numnodes_C;

  value = PLA_Temp_comm_all_info( (PLA_Template) *templ, 
				  &comm_C,
				  &rank_C, &numnodes_C );

  *rank =     ( FINT ) rank_C;
  *numnodes = ( FINT ) numnodes_C;

  *comm = MPI_Comm_c2f( comm_C );

  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMALL_F
#else
void PLA_FOR2C( pla_tempcommall_f )
#endif
( FINT *templ, MPI_Fint *comm, FINT *ierror )
{
  MPI_Comm 
    comm_C = MPI_COMM_NULL;
  int 
    value;

  value = PLA_Temp_comm_all( (PLA_Template) *templ, 
		             &comm_C );

  *comm = MPI_Comm_c2f( comm_C );

  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMALLRANK_F
#else
void PLA_FOR2C( pla_tempcommallrank_f )
#endif
( FINT *templ, FINT *rank, FINT *ierror )
{
  int
    rank_C;

  *ierror = ( FINT ) PLA_Temp_comm_all_rank( (PLA_Template) *templ, 
				             &rank_C );

  *rank = ( FINT ) rank_C;

  return;
}


#if MANUFACTURE == CRAY
void PLA_TEMPCOMMALLSIZE_F
#else
void PLA_FOR2C( pla_tempcommallsize_f )
#endif
( FINT *templ, FINT *numnodes, FINT *ierror )
{
  int
    numnodes_C;

  *ierror = ( FINT ) PLA_Temp_comm_all_size( (PLA_Template) *templ, 
				             &numnodes_C );

  *numnodes = ( FINT ) numnodes_C;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMROWINFO_F
#else
void PLA_FOR2C( pla_tempcommrowinfo_f )
#endif
( FINT *templ, MPI_Fint *comm, FINT *rank, 
			       FINT *numnodes, FINT *ierror )
{
  MPI_Comm 
    comm_C = MPI_COMM_NULL;
  int 
    value,
    rank_C, numnodes_C;

  value = PLA_Temp_comm_row_info( (PLA_Template) *templ, 
				  &comm_C,
				  &rank_C, &numnodes_C );

  *rank = ( FINT ) rank_C;
  *numnodes = ( FINT ) numnodes_C;

  *comm = MPI_Comm_c2f( comm_C );

  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMROW_F
#else
void PLA_FOR2C( pla_tempcommrow_f )
#endif
( FINT *templ, MPI_Fint *comm, FINT *ierror )
{
  MPI_Comm 
    comm_C = MPI_COMM_NULL;
  int 
    value;

  value = PLA_Temp_comm_row( (PLA_Template) *templ, 
		             &comm_C );

  *comm = MPI_Comm_c2f( comm_C );

  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMROWRANK_F
#else
void PLA_FOR2C( pla_tempcommrowrank_f )
#endif
( FINT *templ, FINT *rank, FINT *ierror )
{
  int
    rank_C;

  *ierror = ( FINT ) PLA_Temp_comm_row_rank( ( PLA_Template) *templ, 
					     &rank_C );
  
  *rank = ( FINT ) rank_C;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMROWSIZE_F
#else
void PLA_FOR2C( pla_tempcommrowsize_f )
#endif
( FINT *templ, FINT *numnodes, FINT *ierror )
{
  int
    numnodes_C;

  *ierror = ( FINT ) PLA_Temp_comm_row_size( (PLA_Template) *templ, 
					     &numnodes_C );

  *numnodes = ( FINT ) numnodes_C;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMCOLINFO_F
#else
void PLA_FOR2C( pla_tempcommcolinfo_f )
#endif
( FINT *templ, MPI_Fint *comm, FINT *rank, 
			       FINT *numnodes, FINT *ierror )
{
  MPI_Comm 
    comm_C = MPI_COMM_NULL;
  int 
    value,
    rank_C, numnodes_C;

  value = PLA_Temp_comm_col_info( (PLA_Template) *templ, 
				  &comm_C,
				  &rank_C, &numnodes_C );

  *rank = ( FINT ) rank_C;
  *numnodes = ( FINT ) numnodes_C;

  *comm = MPI_Comm_c2f( comm_C );

  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMCOL_F
#else
void PLA_FOR2C( pla_tempcommcol_f )
#endif
( FINT *templ, MPI_Fint *comm, FINT *ierror )
{
  MPI_Comm 
    comm_C = MPI_COMM_NULL;
  int 
    value;

  value = PLA_Temp_comm_col( (PLA_Template) *templ, 
		             &comm_C );

  *comm = MPI_Comm_c2f( comm_C );

  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMCOLRANK_F
#else
void PLA_FOR2C( pla_tempcommcolrank_f )
#endif
( FINT *templ, FINT *rank, FINT *ierror )
{
  int 
    rank_C;

  *ierror = ( FINT ) PLA_Temp_comm_col_rank( ( PLA_Template ) *templ, 
					     &rank_C );

  *rank = ( FINT ) rank_C;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPCOMMCOLSIZE_F
#else
void PLA_FOR2C( pla_tempcommcolsize_f )
#endif
( FINT *templ, FINT *numnodes, FINT *ierror )
{
  int 
    numnodes_C;

  *ierror = ( FINT ) PLA_Temp_comm_col_size( ( PLA_Template ) *templ, 
				             &numnodes_C );

  *numnodes = ( FINT ) numnodes_C;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPNB_F
#else
void PLA_FOR2C( pla_tempnb_f )
#endif
( FINT *templ, FINT *nb_distr, FINT *ierror )
{
  int
    nb_distr_C;

  *ierror = ( FINT ) PLA_Temp_nb( ( PLA_Template ) *templ, 
		         &nb_distr_C );
  
  *nb_distr = ( FINT ) nb_distr_C;

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPZEROORONE_F
#else
void PLA_FOR2C( pla_tempzeroorone_f )
#endif
( FINT *templ, FINT *zero_or_one, FINT *ierror )
{
  int
    zero_or_one_C;

  *ierror = ( FINT ) PLA_Temp_zero_or_one( ( PLA_Template ) *templ, 
					   &zero_or_one_C );

  *zero_or_one = ( FINT ) zero_or_one_C;

  return;
}

/*
   Section 2.3: Linear Algebra Objects
*/

#if MANUFACTURE == CRAY
void PLA_VECTORCREATE_F
#else
void PLA_FOR2C( pla_vectorcreate_f )
#endif
( FINT *datatype, FINT *m, FINT *templ, 
			  FINT *align_row, FINT *x, FINT *ierror )
{
  printf("Please use pla_mvector_create_f instead\n");
  exit( 0 );

  return;
}

#if MANUFACTURE == CRAY
void PLA_MVECTORCREATE_F
#else
void PLA_FOR2C( pla_mvectorcreate_f )
#endif
( FINT *datatype, FINT *m, FINT *n, FINT *templ, 
                        FINT *align_row, FINT *x, FINT *ierror )
{
  MPI_Datatype 
    datatype_C;
  
  datatype_C = pla_type_f2c( *datatype );

  *ierror = ( FINT ) PLA_Mvector_create( datatype_C, ( int ) *m, ( int ) *n, 
					 ( PLA_Template ) *templ, 
					 ( int ) *align_row, ( PLA_Obj * ) x );

  return;
}

#if MANUFACTURE == CRAY
void PLA_MATRIXCREATE_F
#else
void PLA_FOR2C( pla_matrixcreate_f )
#endif
( FINT *datatype, FINT *m, FINT *n, FINT *templ, 
			  FINT *align_row, FINT *align_col, FINT *a, 
			  FINT *ierror )
{
  MPI_Datatype 
    datatype_C;
  
  datatype_C = pla_type_f2c( *datatype );

  *ierror = ( FINT ) PLA_Matrix_create( datatype_C, ( int ) *m, ( int ) *n, 
					( PLA_Template ) *templ, 
					( int ) *align_row, 
					( int ) *align_col, 
				 	( PLA_Obj * ) a );

  return;
}

#if MANUFACTURE == CRAY
void PLA_MSCALARCREATE_F
#else
void PLA_FOR2C( pla_mscalarcreate_f )
#endif
( FINT *datatype, FINT *owner_row, FINT *owner_col,
			   FINT *m, FINT *n, FINT *templ, FINT *x, 
			   FINT *ierror )
{
  MPI_Datatype 
    datatype_C;
  
  datatype_C = pla_type_f2c( *datatype );

  *ierror = ( FINT ) PLA_Mscalar_create( datatype_C, 
                             ( int ) *owner_row, ( int ) *owner_col,
                             ( int ) *m, ( int ) *n, 
                             ( PLA_Template ) *templ, 
                             ( PLA_Obj * ) x );

  return;
}

#if MANUFACTURE == CRAY
void PLA_PVECTORCREATE_F
#else
void PLA_FOR2C( pla_pvectorcreate_f )
#endif
( FINT *datatype, FINT *project_onto,
                          FINT *owner, FINT *m, FINT *templ, 
                          FINT *align, FINT *x, FINT *ierror )
{
  printf("Please use pla_pmvector_create_f instead\n");
  exit( 0 );
}

#if MANUFACTURE == CRAY
void PLA_PMVECTORCREATE_F
#else
void PLA_FOR2C( pla_pmvectorcreate_f )
#endif
( FINT *datatype, FINT *project_onto,
			    FINT *owner, FINT *m, FINT *n, FINT *templ, 
			    FINT *align, FINT *x, FINT *ierror )
{
  MPI_Datatype 
    datatype_C;
  
  datatype_C = pla_type_f2c( *datatype );

  *ierror = ( FINT ) PLA_Pmvector_create( datatype_C, 
			     ( int ) *project_onto, ( int ) *owner,
                             ( int ) *m, ( int ) *n, 
                             ( PLA_Template ) *templ, 
                             ( int ) *align, (PLA_Obj *) x );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJOBJTYPE_F
#else
void PLA_FOR2C( pla_objobjtype_f )
#endif
( FINT *obj, FINT *objtype, FINT *ierror )
{
  int 
    value;

  *ierror = ( FINT ) PLA_Obj_objtype( (PLA_Obj) *obj, &value );
  *objtype = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJDATATYPE_F
#else
void PLA_FOR2C( pla_objdatatype_f )
#endif
( FINT *obj, FINT *datatype, FINT *ierror )
{
  MPI_Datatype 
    datatype_C;
  int 
    value;
  
  value = PLA_Obj_datatype( (PLA_Obj) *obj, &datatype_C );

  *datatype = pla_type_c2f( datatype_C );
  
  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJTEMPLATE_F
#else
void PLA_FOR2C( pla_objtemplate_f )
#endif
( FINT *obj, FINT *templ, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_template( (PLA_Obj) *obj, (PLA_Template *) templ );

  return;
}

/*
   Section 2.3.2 
*/

#if MANUFACTURE == CRAY
void PLA_OBJFREE_F
#else
void PLA_FOR2C( pla_objfree_f )
#endif
( FINT *obj, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_free( (PLA_Obj *) obj );

  return;
}

/*
   Section 2.3.3 
*/

#if MANUFACTURE == CRAY
void PLA_OBJGLOBALINFO_F
#else
void PLA_FOR2C( pla_objglobalinfo_f )
#endif
( FINT *obj, FINT *global_length,
			   FINT *global_width, FINT *project_onto,
			   FINT *owner_row, FINT *owner_col,
			   FINT *global_align_row, FINT *global_align_col, FINT *ierror )
{
  int
    global_length_C, global_width_C, project_onto_C, owner_row_C, 
    owner_col_C, global_align_row_C, global_align_col_C;

  *ierror = ( FINT ) PLA_Obj_global_info( (PLA_Obj) *obj, &global_length_C,
			   &global_width_C, &project_onto_C,
			   &owner_row_C, &owner_col_C,
			   &global_align_row_C, &global_align_col_C );

  *global_length = ( FINT ) global_length_C; 
  *global_width = ( FINT ) global_width_C; 
  *project_onto = ( FINT ) project_onto_C; 
  *owner_row = ( FINT ) owner_row_C; 
  *owner_col = ( FINT ) owner_col_C; 
  *global_align_row = ( FINT ) global_align_row_C; 
  *global_align_col = ( FINT ) global_align_col_C;
  
  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJGLOBALLENGTH_F
#else
void PLA_FOR2C( pla_objgloballength_f )
#endif
( int *obj, int *length, int *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_global_length( (PLA_Obj) *obj, &value );

  *length = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJGLOBALWIDTH_F
#else
void PLA_FOR2C( pla_objglobalwidth_f )
#endif
( FINT *obj, FINT *width, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_global_width( (PLA_Obj) *obj, &value );

  *width = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJGLOBALSIZE_F
#else
void PLA_FOR2C( pla_objglobalsize_f )
#endif
( FINT *obj, FINT *size, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_global_size( (PLA_Obj) *obj, &value );

  *size = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJGLOBALNUMVECS_F
#else
void PLA_FOR2C( pla_objglobalnumvecs_f )
#endif
( FINT *obj, FINT *numvecs, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_global_numvecs( (PLA_Obj) *obj, &value );

  *numvecs = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJPROJECTONTO_F
#else
void PLA_FOR2C( pla_objprojectonto_f )
#endif
( FINT *obj, FINT *project_onto, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_project_onto( (PLA_Obj) *obj, &value );

  *project_onto = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJOWNERROW_F
#else
void PLA_FOR2C( pla_objownerrow_f )
#endif
( FINT *obj, FINT *owner_row, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_owner_row( (PLA_Obj) *obj, &value );

  *owner_row = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJOWNERCOL_F
#else
void PLA_FOR2C( pla_objownercol_f )
#endif
( FINT *obj, FINT *owner_col, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_owner_col( (PLA_Obj) *obj, &value );

  *owner_col = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJGLOBALALIGN_F
#else
void PLA_FOR2C( pla_objglobalalign_f )
#endif
( FINT *obj, FINT *global_align, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_global_align( (PLA_Obj) *obj, &value );

  *global_align = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJGLOBALALIGNROW_F
#else
void PLA_FOR2C( pla_objglobalalignrow_f )
#endif
( FINT *obj, FINT *global_align_row, 
				 FINT *ierror )
{
  int
    value;
 
  *ierror = ( FINT ) PLA_Obj_global_align_row( (PLA_Obj) *obj, &value );

  *global_align_row = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJGLOBALALIGNCOL_F
#else
void PLA_FOR2C( pla_objglobalaligncol_f )
#endif
( FINT *obj, FINT *global_align_col, 
				 FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_global_align_col( (PLA_Obj) *obj, &value );

  *global_align_col = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJLOCALINFO_F
#else
void PLA_FOR2C( pla_objlocalinfo_f )
#endif
( FINT *obj, FINT *local_length,
			  FINT *local_width, void **dummy,
			  FINT *local_stride, FINT *local_ldim, FINT *ierror )
{
  int 
    value,
    local_length_C,
    local_width_C, 
    local_stride_C, 
    local_ldim_C;

  value = PLA_Obj_local_info( (PLA_Obj) *obj, &local_length_C,
			  &local_width_C, dummy,
			  &local_stride_C, &local_ldim_C );

  *dummy = NULL;
  *local_length = ( FINT ) local_length_C;
  *local_width = ( FINT ) local_width_C;
  *local_stride = ( FINT ) local_stride_C;
  *local_ldim = ( FINT ) local_ldim_C;

  *ierror = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJLOCALLENGTH_F
#else
void PLA_FOR2C( pla_objlocallength_f )
#endif
( FINT *obj, FINT *local_length, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_local_length( (PLA_Obj) *obj, &value );

  *local_length = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJLOCALWIDTH_F
#else
void PLA_FOR2C( pla_objlocalwidth_f )
#endif
( FINT *obj, FINT *local_width, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_local_width( (PLA_Obj) *obj, &value );

  *local_width = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJLOCALBUFFER_F
#else
void PLA_FOR2C( pla_objlocalbuffer_f )
#endif
( FINT *obj, void **local_buffer, FINT *ierror )
{
  printf("pla_obj_local_buffer cannot be called from FORTRAN\n");
  exit( 0 );
}

#if MANUFACTURE == CRAY
void PLA_OBJLOCALSTRIDE_F
#else
void PLA_FOR2C( pla_objlocalstride_f )
#endif
( FINT *obj, FINT *local_stride, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_local_stride( (PLA_Obj) *obj, &value );

  *local_stride = ( FINT ) value;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJLOCALLDIM_F
#else
void PLA_FOR2C( pla_objlocalldim_f )
#endif
( FINT *obj, FINT *local_ldim, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_local_ldim( (PLA_Obj) *obj, &value );

  *local_ldim = ( FINT ) value;

  return;
}

/*
   Section 2.3.4
*/

#if MANUFACTURE == CRAY
void PLA_OBJGETLOCALCONTENTS_F
#else
void PLA_FOR2C( pla_objgetlocalcontents_f )
#endif
( FINT *obj, FINT *trans, 
				  FINT *m, FINT *n, void *address,
				  FINT *lda, FINT *stride, FINT *ierror )
{
  int 
    m_C, n_C;

  *ierror = ( FINT ) PLA_Obj_get_local_contents( ( PLA_Obj ) *obj, 
						 ( int ) *trans, 
	                                         &m_C, &n_C, 
						 address, ( int ) *lda, 
						 ( int ) *stride );
  *m = ( FINT ) m_C;
  *n = ( FINT ) n_C;

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJSETLOCALCONTENTS_F
#else
void PLA_FOR2C( pla_objsetlocalcontents_f )
#endif
( FINT *trans, 
				  FINT *m, FINT *n, void *address,
				  FINT *lda, FINT *stride,
				  FINT *obj, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_set_local_contents( ( int ) *trans, 
						 ( int ) *m, 
						 ( int ) *n, 
						 address, 
						 ( int ) *lda, 
						 ( int ) *stride,
						 (PLA_Obj) *obj );

  return;
}

/*
   Section 2.3.5
*/

#if MANUFACTURE == CRAY
void PLA_OBJSET_F
#else
void PLA_FOR2C( pla_objset_f )
#endif
( FINT *obj, FINT *datatype, void *value, FINT *ierror )
{
  MPI_Datatype datatype_C;
  
  datatype_C = pla_type_f2c( *datatype );

  *ierror = ( FINT ) PLA_Obj_set( (PLA_Obj) *obj, datatype_C, value );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJSETTOONE_F
#else
void PLA_FOR2C( pla_objsettoone_f )
#endif
( FINT *obj, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_set_to_one( (PLA_Obj) *obj );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJSETTOZERO_F
#else
void PLA_FOR2C( pla_objsettozero_f )
#endif
( int *obj, int *ierror )
{
  *ierror = ( FINT ) PLA_Obj_set_to_zero( (PLA_Obj) *obj );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJSETTOMINUSONE_F
#else
void PLA_FOR2C( pla_objsettominusone_f )
#endif
( int *obj, int *ierror )
{
  *ierror = ( FINT ) PLA_Obj_set_to_minus_one( (PLA_Obj) *obj );

  return;
}

/*************************************************************
   Chapter 3
**************************************************************/

/* 
   Section 3.1: Creating Views into Objects
*/

#if MANUFACTURE == CRAY
void PLA_OBJVIEW_F
#else
void PLA_FOR2C( pla_objview_f )
#endif
( FINT *obj, FINT *global_length, 
		    FINT *global_width, FINT *align_row, 
		    FINT *align_col, FINT *obj_new, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_view( ( PLA_Obj ) *obj, 
				   ( int ) *global_length, 
				   ( int ) *global_width, 
				   ( int ) *align_row, 
				   ( int ) *align_col, 
				   (PLA_Obj *) obj_new );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJVIEWALL_F
#else
void PLA_FOR2C( pla_objviewall_f )
#endif
( FINT *obj, FINT *obj_new, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_view_all( (PLA_Obj) *obj, (PLA_Obj *) obj_new );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJVIEWSWAP_F
#else
void PLA_FOR2C( pla_objviewswap_f )
#endif
( FINT *obj1, FINT *obj2, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_view_swap( (PLA_Obj *) obj1, (PLA_Obj *) obj2 );

  return;
}

/* 
   Section 3.2: Splitting of Linear Algebra Objects
*/

#if MANUFACTURE == CRAY
void PLA_OBJSPLIT4_F
#else
void PLA_FOR2C( pla_objsplit4_f )
#endif
( FINT *obj, FINT *size_row, FINT *size_col, 
			FINT *obj_11, FINT *obj_12,
			FINT *obj_21, FINT *obj_22, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_split_4( (PLA_Obj) *obj, 
				      ( int ) *size_row, ( int ) *size_col, 
				      (PLA_Obj *) obj_11, (PLA_Obj *) obj_12, 
				      (PLA_Obj *) obj_21, (PLA_Obj *) obj_22 );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJHORZSPLIT2_F
#else
void PLA_FOR2C( pla_objhorzsplit_2_f )
#endif
( FINT *obj, FINT *size, FINT *obj_1, FINT *obj_2, 
			     FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_horz_split_2( (PLA_Obj) *obj, 
					   ( int ) *size, 
					   (PLA_Obj *) obj_1, 
					   (PLA_Obj *) obj_2 ); 

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJVERTSPLIT2_F
#else
void PLA_FOR2C( pla_obj_vert_split_2_f )
#endif
( FINT *obj, FINT *size, FINT *obj_1, FINT *obj_2, 
			     FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_vert_split_2( (PLA_Obj) *obj, 
					   ( int ) *size, 
					   (PLA_Obj *) obj_1, 
					   (PLA_Obj *) obj_2 ); 

  return;
}

/* 
   Section 3.3: Shifting of Linear Algebra Objects
*/

#if MANUFACTURE == CRnAY
void PLA_OBJVIEWSHIFT_F
#else
void PLA_FOR2C( pla_objviewshift_f )
#endif
( FINT *obj,       FINT *top, 
                                FINT *left,         FINT *right,
                                          FINT *bottom, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_view_shift( ( PLA_Obj ) *obj,       
					 ( int ) *top, ( int ) *left,  
					 ( int ) *right, ( int ) *bottom );

  return;
}

/* 
   Section 3.4: Determining Where to Split
*/

#if MANUFACTURE == CRAY
void PLA_OBJSPLITSIZE_F
#else
void PLA_FOR2C( pla_objsplitsize_f )
#endif
( FINT *obj, FINT *side, FINT *size, FINT *owner, 
			   FINT *ierror )
{
  int 
    size_C, owner_C;

  *ierror = ( FINT ) PLA_Obj_split_size( (PLA_Obj) *obj, 
					 ( int ) *side, &size_C, &owner_C );

  *size = ( FINT ) size_C;
  *owner = ( FINT ) owner_C;

  return;
}

/* 
   Section 3.5: Creating Objects "Conformal to" Other Objects
*/

#if MANUFACTURE == CRAY
void PLA_VECTORCREATECONFTO_F
#else
void PLA_FOR2C( pla_vectorcreateconfnto_f )
#endif
( FINT *obj,  FINT *new_obj, FINT *ierror )
{
  printf("Please use pla_mvector_create_conf_to_f instead of\n");
  printf("           pla_vector_create_conf_to_f\n");
  exit( 0 );
}

#if MANUFACTURE == CRAY
void PLA_MVECTORCREATECONFTO_F
#else
void PLA_FOR2C( pla_mvectorcreateconfto_f )
#endif
( FINT *obj, FINT *num_vectors, FINT *new_obj,
				   FINT *ierror )
{
  *ierror = ( FINT ) PLA_Mvector_create_conf_to( (PLA_Obj) *obj, 
						 ( int ) *num_vectors, 
						 (PLA_Obj *) new_obj );

  return;
}

#if MANUFACTURE == CRAY
void PLA_PVECTORCREATECONFTO_F
#else
void PLA_FOR2C( pla_pvectorcreateconfto_f )
#endif
( FINT *obj,  FINT *project_onto,
				   FINT *owner, FINT *new_obj, FINT *ierror )
{
  printf("Please use pla_pmvector_create_conf_to_f instead of\n");
  printf("           pla_pvector_create_conf_to_f\n");
  exit( 0 );
}

#if MANUFACTURE == CRAY
void PLA_PMVECTORCREATECONFTO_F
#else
void PLA_FOR2C( pla_pmvectorcreateconfto_f )
#endif
( FINT *obj, FINT *project_onto,
				    FINT *owner, FINT *num_vectors, 
				    FINT *new_obj, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Pmvector_create_conf_to( (PLA_Obj) *obj, 
						  ( int ) *project_onto, 
						  ( int ) *owner, 
						  ( int ) *num_vectors, 
						  ( PLA_Obj * ) new_obj );

  return;
}

#if MANUFACTURE == CRAY
void PLA_MATRIXCREATECONFTO_F
#else
void PLA_FOR2C( pla_matrixcreateconfto_f )
#endif
( FINT *obj,  FINT *new_obj, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Matrix_create_conf_to( ( PLA_Obj ) *obj,  
						( PLA_Obj * ) new_obj );

  return;
}

#if MANUFACTURE == CRAY
void PLA_MSCALARCREATECONFTO_F
#else
void PLA_FOR2C( pla_mscalarcreateconfto_f )
#endif
( FINT *obj,  FINT *owner_row, 
				   FINT *owner_col, FINT *new_obj, 
				   FINT *ierror )
{
  *ierror = ( FINT ) PLA_Mscalar_create_conf_to( ( PLA_Obj ) *obj, 
						 ( int ) *owner_row,
						 ( int ) *owner_col, 
						 ( PLA_Obj * ) new_obj );

  return;
}

/* 
   Section 3.6: Annotating Object Orientation
*/

#if MANUFACTURE == CRAY
void PLA_OBJSETORIENTATION_F
#else
void PLA_FOR2C( pla_objsetorientation_f )
#endif
( FINT *obj,  FINT *project_onto, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_set_orientation( ( PLA_Obj ) *obj, 
					      ( int ) *project_onto );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJGETORIENTATION_F
#else
void PLA_FOR2C( pla_objgetorientation_f )
#endif
( FINT *obj,  FINT *project_onto, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_get_orientation( (PLA_Obj) *obj, &value );

  *project_onto = ( FINT ) value;

  return;
}

/* 
   Section 3.7: Casting Object Types
*/

#if MANUFACTURE == CRAY
void PLA_OBJOBJTYPECAST_F
#else
void PLA_FOR2C( pla_objobjtypecast_f )
#endif
( FINT *obj,  FINT *objtype, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_objtype_cast( (PLA_Obj) *obj, ( int ) *objtype );

  return;
}

/*************************************************************
   Chapter 4
**************************************************************/

/* 
   Section 4.1: Introduction
*/

/* 
   Section 4.2: API-Activation
*/

#if MANUFACTURE == CRAY
void PLA_APIBEGIN_F
#else
void PLA_FOR2C( pla_apibegin_f )
#endif
( FINT *ierror )
{
  *ierror = ( FINT ) PLA_API_begin( );

  return;
}

#if MANUFACTURE == CRAY
void PLA_APIEND_F
#else
void PLA_FOR2C( pla_apiend_f )
#endif
( FINT *ierror )
{
  *ierror = ( FINT ) PLA_API_end( );

  return;
}

#if MANUFACTURE == CRAY
void PLA_APISTATE_F
#else
void PLA_FOR2C( pla_apistate_f )
#endif
( FINT *state, FINT *ierror )
{
  int 
    value;
  
  *ierror = ( FINT ) PLA_API_state( &value );
  
  *state = ( FINT ) value;

  return;
}

/* 
   Section 4.3: Opening and Closing Objects
*/

#if MANUFACTURE == CRAY
void PLA_OBJAPIOPEN_F
#else
void PLA_FOR2C( pla_objapiopen_f )
#endif
( FINT *obj, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_API_open( ( PLA_Obj ) *obj );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJAPICLOSE_F
#else
void PLA_FOR2C( pla_objapiclose_f )
#endif
( FINT * obj, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_API_close( (PLA_Obj) *obj );

  return;
}

#if MANUFACTURE == CRAY
void PLA_OBJAPIMODEF
#else
void PLA_FOR2C( pla_objapimodef )
#endif
( FINT *obj, FINT *mode, FINT *ierror )
{
  int
    value;

  *ierror = ( FINT ) PLA_Obj_API_mode( (PLA_Obj) *obj, &value );

  *mode = ( FINT ) value;

  return;
}

/* 
   Section 4.4: Accessing a Vector
*/

#if MANUFACTURE == CRAY
void PLA_APIAXPYVECTORTOGLOBAL_F
#else
void PLA_FOR2C( pla_apiaxpyvectortoglobal_f )
#endif
( FINT *size, 
				     void *alpha, void *x, 
				     FINT *stride, FINT *x_global, 
				     FINT *disp_row, FINT *ierror )
{
  *ierror = ( FINT ) PLA_API_axpy_vector_to_global( ( int ) *size,
				   alpha, x, 
				   ( int ) *stride, (PLA_Obj) *x_global, 
				   ( int ) *disp_row );

  return;
}

#if MANUFACTURE == CRAY
void PLA_APIAXPYGLOBALTOVECTOR_F
#else
void PLA_FOR2C( pla_apiaxpyglobaltovector_f )
#endif
( FINT *size, 
				     void *alpha, FINT *x_global,
				     FINT *disp_row, 
                                     void *x, FINT *stride, FINT *ierror )
{
  *ierror = ( FINT ) PLA_API_axpy_global_to_vector( ( int ) *size,
				   alpha, (PLA_Obj) *x_global, 
				   ( int ) *disp_row, 
				   x, ( int ) *stride ); 

  return;
}

#if MANUFACTURE == CRAY
void PLA_APIMULTIAXPYVECTORTOGLOBAL_F
#else
void PLA_FOR2C( pla_apimultiaxpyvectortoglobal_f )
#endif
(  FINT *nsub, FINT *sizes, 
				     void *alpha, void *local_vector, 
				     FINT *stride, FINT *x_global, 
				     FINT *disps_row, FINT *ierror )
{
  int 
    *sizes_C, *disps_row_C, i;

  sizes_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub ) * sizeof( int ) );
  disps_row_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub ) * sizeof( int ) );

  for ( i=0; i< ( int ) *nsub; i++ ){
    sizes_C[ i ] = ( int ) sizes[ i ];
    disps_row_C[ i ] = ( int ) disps_row[ i ];
  }

  *ierror = ( FINT ) PLA_API_multi_axpy_vector_to_global(
                                   ( int ) *nsub, sizes_C,
				   alpha, local_vector, 
				   ( int ) *stride, (PLA_Obj) *x_global, 
				   disps_row_C );

  PLA_free( sizes_C );
  PLA_free( disps_row_C );

  return;
}

#if MANUFACTURE == CRAY
void PLA_APIMULTIAXPYGLOBALTOVECTOR_F
#else
void PLA_FOR2C( pla_apimultiaxpyglobaltovector_f )
#endif
( FINT *nsub, FINT *sizes, 
				     void *alpha, FINT *x_global,
				     FINT *disps_row, 
                                     void *x, FINT *stride, FINT *ierror )
{
  int 
    *sizes_C, *disps_row_C, i;

  sizes_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub ) * sizeof( int ) );
  disps_row_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub ) * sizeof( int ) );

  for ( i=0; i< ( int ) *nsub; i++ ){
    sizes_C[ i ] = ( int ) sizes[ i ];
    disps_row_C[ i ] = ( int ) disps_row[ i ];
  }

  *ierror = ( FINT ) PLA_API_multi_axpy_global_to_vector(
				   ( int ) *nsub, 
				   sizes_C,
				   alpha, 
				   (PLA_Obj) *x_global, 
				   disps_row_C, 
				   x, 
				   ( int ) *stride ); 

  PLA_free( sizes_C );
  PLA_free( disps_row_C );

  return;
}

/* 
   Section 4.5:  Accessing a Matrix
*/

#if MANUFACTURE == CRAY
void PLA_APIAXPYMATRIXTOGLOBAL_F
#else
void PLA_FOR2C( pla_apiaxpymatrixtoglobal_f )
#endif
( FINT *size_row, FINT *size_col,
				      void *alpha, void *a, 
				      FINT *lda, FINT *a_global, 
				      FINT *disp_row, FINT *disp_col, 
				      FINT *ierror )
{
  *ierror = ( FINT ) PLA_API_axpy_matrix_to_global(
				   ( int ) *size_row, ( int ) *size_col,
				   alpha, a, 
				   ( int ) *lda, (PLA_Obj) *a_global, 
				   ( int ) *disp_row, ( int ) *disp_col );

  return;
}

#if MANUFACTURE == CRAY
void PLA_APIAXPYGLOBALTOMATRIX_F
#else
void PLA_FOR2C( pla_apiaxpyglobaltomatrix_f )
#endif
( FINT *size_row, FINT *size_col,
				     void *alpha, FINT *a_global, 
				     FINT *disp_row, FINT *disp_col,
				     void *a, FINT *lda, FINT *ierror )
{
  *ierror = ( FINT ) PLA_API_axpy_global_to_matrix( 
				   ( int ) *size_row, ( int ) *size_col,
				   alpha, (PLA_Obj) *a_global, 
				   ( int ) *disp_row, ( int ) *disp_col,
				   a, ( int ) *lda );

  return;
}

#if MANUFACTURE == CRAY
void PLA_APIMULTIAXPYMATRIXTOGLOBAL_F
#else
void PLA_FOR2C( pla_apimultiaxpymatrixtoglobal_f )
#endif
( FINT *nsub_row, FINT *nsub_col,
					    FINT *size_row, FINT *size_col,
					    void *alpha, void *a, 
					    FINT *lda, FINT *a_global, 
					    FINT *disp_row, FINT *disp_col, 
					    FINT *ierror )
{
  int 
    *size_row_C, *disps_row_C,
    *size_col_C, *disps_col_C,
    i;

  size_row_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub_row ) * sizeof( int ) );
  disps_row_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub_row ) * sizeof( int ) );

  size_col_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub_col ) * sizeof( int ) );
  disps_col_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub_col ) * sizeof( int ) );

  for ( i=0; i< ( int ) *nsub_row; i++ ){
    size_row_C[ i ] = ( int ) size_row[ i ];
    disps_row_C[ i ] = ( int ) disp_row[ i ];
  }

  for ( i=0; i< ( int ) *nsub_col; i++ ){
    size_col_C[ i ] = ( int ) size_col[ i ];
    disps_col_C[ i ] = ( int ) disp_col[ i ];
  }

  *ierror = ( FINT ) PLA_API_multi_axpy_matrix_to_global( 
				   ( int ) *nsub_row, ( int ) *nsub_col,
				   size_row_C, size_col_C,
				   alpha, a, 
				   ( int ) *lda, (PLA_Obj) *a_global, 
				   disps_row_C, disps_col_C );

  PLA_free( size_row_C );
  PLA_free( disps_row_C );
  PLA_free( size_col_C );
  PLA_free( disps_col_C );

  return;
}

#if MANUFACTURE == CRAY
void PLA_APIMULTIAXPYGLOBALTOMATRIX_F
#else
void PLA_FOR2C( pla_apimultiaxpyglobaltomatrix_f )
#endif
( FINT *nsub_row, FINT *nsub_col,
			             FINT *size_row, FINT *size_col,
				     void *alpha, FINT *a_global, 
				     FINT *disp_row, FINT *disp_col,
				     void *a, FINT *lda, FINT *ierror )
{
  int 
    *size_row_C, *disps_row_C,
    *size_col_C, *disps_col_C,
    i;

  size_row_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub_row ) * sizeof( int ) );
  disps_row_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub_row ) * sizeof( int ) );

  size_col_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub_col ) * sizeof( int ) );
  disps_col_C = ( int * ) PLA_malloc( (size_t) ( ( int ) *nsub_col ) * sizeof( int ) );

  for ( i=0; i< ( int ) *nsub_row; i++ ){
    size_row_C[ i ] = ( int ) size_row[ i ];
    disps_row_C[ i ] = ( int ) disp_row[ i ];
  }

  for ( i=0; i< ( int ) *nsub_col; i++ ){
    size_col_C[ i ] = ( int ) size_col[ i ];
    disps_col_C[ i ] = ( int ) disp_col[ i ];
  }

  *ierror = ( FINT ) PLA_API_multi_axpy_global_to_matrix( *nsub_row, *nsub_col,
				   size_row, size_col,
				   alpha, (PLA_Obj) *a_global, 
				   disps_row_C, disps_col_C,
				   a, *lda );

  PLA_free( size_row_C );
  PLA_free( disps_row_C );
  PLA_free( size_col_C );
  PLA_free( disps_col_C );

  return;
}

/* 
   Section 4.6: Completion and Synchronization
*/

#if MANUFACTURE == CRAY
void PLA_OBJAPISYNC_F
#else
void PLA_FOR2C( pla_objapisync_f )
#endif
( FINT *obj, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_API_sync( (PLA_Obj) *obj );

  return;
}

/*************************************************************
   Chapter 5
**************************************************************/

/* 
   Section 5.1: Copy
*/

#if MANUFACTURE == CRAY
void PLA_COPY_F
#else
void PLA_FOR2C( pla_copy_f )
#endif
( FINT *obj_from, FINT *obj_to, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Copy( (PLA_Obj) *obj_from, (PLA_Obj) *obj_to );

  return;
}

#if MANUFACTURE == CRAY
void PLA_COPYX_F
#else
void PLA_FOR2C( pla_copyx_f )
#endif
( FINT *shape, FINT *obj_from, FINT *obj_to, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Copy_x( ( int ) *shape, (PLA_Obj) *obj_from, 
				 (PLA_Obj) *obj_to ); 

  return;
}

/* 
   Section 5.2:
*/

#if MANUFACTURE == CRAY
void PLA_REDUCE_F
#else
void PLA_FOR2C( pla_reduce_f )
#endif
( FINT *obj_from, FINT *op, FINT *obj_to, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Reduce( (PLA_Obj) *obj_from, MPI_Op_f2c( *op), 
		     (PLA_Obj) *obj_to );

  return;
}


#if MANUFACTURE == CRAY
void PLA_REDUCEX_F
#else
void PLA_FOR2C( pla_reducex_f )
#endif
( FINT *shape, FINT *obj_from, FINT *obj_alpha, 
		    FINT *obj_to, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Reduce_x( ( int ) *shape, (PLA_Obj) *obj_from, 
				   (PLA_Obj) *obj_alpha,
				   (PLA_Obj) *obj_to );

  return;
}

/* 
   Section 5.3: Pipelining Computation and Communication
*/

#if MANUFACTURE == CRAY
void PLA_TEMPSETCOMMDIR_F
#else
void PLA_FOR2C( pla_tempsetcommdir_f )
#endif
( FINT *templ, FINT *mesh_dimension, 
			      FINT *direction, FINT *ierror )
{
  printf("pla_temp_set_comm_dir_f not currently available\n");
  exit( 0 );
/*  *ierror = ( FINT ) PLA_Temp_set_comm_dir( (PLA_Template) *templ, *mesh_dimension, 
			        *direction ); */

  return;
}

#if MANUFACTURE == CRAY
void PLA_TEMPGETCOMMDIR_F
#else
void PLA_FOR2C( pla_tempgetcommdir_f )
#endif
( FINT *templ, FINT *mesh_dimension, 
			      FINT *direction, FINT *ierror )
{
  printf("pla_temp_get_comm_dir_f not currently available\n");
  exit( 0 );
/*  *ierror = ( FINT ) PLA_Temp_get_comm_dir( (PLA_Template) *templ, *mesh_dimension, 
			        direction ); */

  return;
}

/* 
   Section 5.4:
*/


/*************************************************************
   Chapter 6
**************************************************************/

/*
   Local level-1 BLAS
*/

#if MANUFACTURE == CRAY
void PLA_LOCALCOPY_F
#else
void PLA_FOR2C( pla_localcopy_f )
#endif
( FINT *x, FINT *y, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_copy( (PLA_Obj) *x, (PLA_Obj) *y );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALSWAP_F
#else
void PLA_FOR2C( pla_localswap_f )
#endif
( FINT *x, FINT *y, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_swap( (PLA_Obj) *x, (PLA_Obj) *y );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALSCAL_F
#else
void PLA_FOR2C( pla_localscal_f )
#endif
( FINT *alpha, FINT *x, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_scal( (PLA_Obj) *alpha, (PLA_Obj) *x );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALINVSCAL_F
#else
void PLA_FOR2C( pla_localinvscal_f )
#endif
( FINT *alpha, FINT *x, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_inv_scal( (PLA_Obj) *alpha, (PLA_Obj) *x );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALAXPY_F
#else
void PLA_FOR2C( pla_localaxpy_f )
#endif
( FINT *alpha, FINT *x, FINT *y, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_axpy( (PLA_Obj) *alpha, 
		   (PLA_Obj) *x, 
		   (PLA_Obj) *y );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALDOT_F
#else
void PLA_FOR2C( pla_localdot_f )
#endif
( FINT *x, FINT *y, FINT *alpha , FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_dot( (PLA_Obj) *x, (PLA_Obj) *y, (PLA_Obj) *alpha );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALASUM_F
#else
void PLA_FOR2C( pla_localasum_f )
#endif
( FINT *x, FINT *alpha , FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_asum( (PLA_Obj) *x, (PLA_Obj) *alpha );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALNRM2_F
#else
void PLA_FOR2C( pla_localnrm2_f )
#endif
( FINT *x, FINT *alpha , FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_nrm2( (PLA_Obj) *x, (PLA_Obj) *alpha );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALIAMAX_F
#else
void PLA_FOR2C( pla_localiamax_f )
#endif
( FINT *x, FINT *alpha, FINT *index, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_iamax( (PLA_Obj) *x,
		   (PLA_Obj) *alpha, 
		   (PLA_Obj) *index );

  return;
}

/*
   Global level-1 BLAS
*/

#if MANUFACTURE == CRAY
void PLA_SWAP_F
#else
void PLA_FOR2C( pla_swap_f )
#endif
( FINT *x, FINT *y, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Swap( (PLA_Obj) *x, (PLA_Obj) *y );

  return;
}

#if MANUFACTURE == CRAY
void PLA_SCAL_F
#else
void PLA_FOR2C( pla_scal_f )
#endif
( FINT *alpha, FINT *x, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Scal( (PLA_Obj) *alpha, (PLA_Obj) *x );

  return;
}

#if MANUFACTURE == CRAY
void PLA_INVSCAL_F
#else
void PLA_FOR2C( pla_invscal_f )
#endif
( FINT *alpha, FINT *x, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Inv_scal( (PLA_Obj) *alpha, (PLA_Obj) *x );

  return;
}

#if MANUFACTURE == CRAY
void PLA_AXPY_F
#else
void PLA_FOR2C( pla_axpy_f )
#endif
( FINT *alpha, FINT *x, FINT *y, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Axpy( (PLA_Obj) *alpha, 
		   (PLA_Obj) *x, 
		   (PLA_Obj) *y );

  return;
}

#if MANUFACTURE == CRAY
void PLA_DOT_F
#else
void PLA_FOR2C( pla_dot_f )
#endif
( FINT *x, FINT *y, FINT *alpha , FINT *ierror )
{
  *ierror = ( FINT ) PLA_Dot( (PLA_Obj) *x, (PLA_Obj) *y, (PLA_Obj) *alpha );

  return;
}

#if MANUFACTURE == CRAY
void PLA_ASUM_F
#else
void PLA_FOR2C( pla_asum_f )
#endif
( FINT *x, FINT *alpha , FINT *ierror )
{
  *ierror = ( FINT ) PLA_Asum( (PLA_Obj) *x, (PLA_Obj) *alpha );

  return;
}

#if MANUFACTURE == CRAY
void PLA_NRM2_F
#else
void PLA_FOR2C( pla_nrm2_f )
#endif
( FINT *x, FINT *alpha , FINT *ierror )
{
  *ierror = ( FINT ) PLA_Nrm2( (PLA_Obj) *x, (PLA_Obj) *alpha );

  return;
}

#if MANUFACTURE == CRAY
void PLA_IAMAX_F
#else
void PLA_FOR2C( pla_iamax_f )
#endif
( FINT *x, FINT *index, FINT *alpha, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Iamax( (PLA_Obj) *x,
		   (PLA_Obj) *index, 
		   (PLA_Obj) *alpha );

  return;
}

/*************************************************************
   Chapter 7
**************************************************************/

/*
   Local level-2 BLAS
*/

#if MANUFACTURE == CRAY
void PLA_LOCALGEMV_F
#else
void PLA_FOR2C( pla_localgemv_f )
#endif
( FINT *trans, FINT *alpha, FINT *a, FINT *x, 
		                    FINT *beta, FINT *y, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_gemv( ( int ) *trans, 
				     (PLA_Obj) *alpha, (PLA_Obj) *a, 
				     (PLA_Obj) *x, 
				     (PLA_Obj) *beta,  (PLA_Obj) *y );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALSYMV_F
#else
void PLA_FOR2C( pla_localsymv_f )
#endif
( FINT *uplo, FINT *alpha, FINT *a, FINT *x, FINT *beta, FINT *y, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_symv( ( int ) *uplo, 
				     (PLA_Obj) *alpha, (PLA_Obj) *a, 
				     (PLA_Obj) *x, 
				     (PLA_Obj) *beta,  (PLA_Obj) *y );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALGER_F
#else
void PLA_FOR2C( pla_localger_f )
#endif
( FINT *alpha, FINT *x, FINT *y, FINT *a, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_ger( (PLA_Obj) *alpha, 
				    (PLA_Obj) *x, (PLA_Obj) *y, 
				    (PLA_Obj) *a );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALSYR_F
#else
void PLA_FOR2C( pla_localsyr_f )
#endif
( FINT *uplo, FINT *alpha, FINT *x, FINT *a, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_syr( ( int ) *uplo, 
				    (PLA_Obj) *alpha, (PLA_Obj) *x, 
				    (PLA_Obj) *a );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALSYR2_F
#else
void PLA_FOR2C( pla_localsyr2_f )
#endif
( FINT *uplo, FINT *alpha, FINT *x, FINT *y, 
		       FINT *a, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_syr2( ( int ) *uplo, 
				     (PLA_Obj) *alpha, (PLA_Obj) *x, 
				     (PLA_Obj) *y, 
				     (PLA_Obj) *a );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALTRMV_F
#else
void PLA_FOR2C( pla_localtrmv_f )
#endif
( FINT *uplo, FINT *trans, FINT *diag, 
		       FINT *a, FINT *x, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_trmv( ( int ) *uplo, ( int ) *trans, 
				     ( int ) *diag, 
				     (PLA_Obj) *a, (PLA_Obj) *x ); 

  return;
}


#if MANUFACTURE == CRAY
void PLA_LOCALTRSV_F
#else
void PLA_FOR2C( pla_localtrsv_f )
#endif
( FINT *uplo, FINT *trans, FINT *diag, FINT *a, FINT *x, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_trsv( ( int ) *uplo, ( int ) *trans, 
				     ( int ) *diag, 
				     (PLA_Obj) *a, (PLA_Obj) *x ); 

  return;
}

/*
   Global level-2 BLAS
*/

#if MANUFACTURE == CRAY
void PLA_GEMV_F
#else
void PLA_FOR2C( pla_gemv_f )
#endif
( FINT *trans, FINT *alpha, FINT *a, FINT *x, 
		 FINT *beta, FINT *y, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Gemv( ( int ) *trans, 
			       (PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *x, 
			       (PLA_Obj) *beta,  (PLA_Obj) *y );

  return;
}

#if MANUFACTURE == CRAY
void PLA_SYMV_F
#else
void PLA_FOR2C( pla_symv_f )
#endif
( FINT *uplo, FINT *alpha, FINT *a, FINT *x, 
		 FINT *beta, FINT *y, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Symv( ( int )  *uplo, 
			       (PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *x, 
			       (PLA_Obj) *beta,  (PLA_Obj) *y );

  return;
}

#if MANUFACTURE == CRAY
void PLA_GER_F
#else
void PLA_FOR2C( pla_ger_f )
#endif
( FINT *alpha, FINT *x, FINT *y, FINT *a, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Ger( (PLA_Obj) *alpha, (PLA_Obj) *x, (PLA_Obj) *y, 
			      (PLA_Obj) *a );

  return;
}

#if MANUFACTURE == CRAY
void PLA_SYR_F
#else
void PLA_FOR2C( pla_syr_f )
#endif
( FINT *uplo, FINT *alpha, FINT *x, FINT *a, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Syr( ( int ) *uplo, 
			      (PLA_Obj) *alpha, (PLA_Obj) *x, (PLA_Obj) *a );

  return;
}

#if MANUFACTURE == CRAY
void PLA_SYR2_F
#else
void PLA_FOR2C( pla_syr2_f )
#endif
( FINT *uplo, FINT *alpha, FINT *x, FINT *y, FINT *a, 
		 FINT *ierror )
{
  *ierror = ( FINT ) PLA_Syr2( ( int ) *uplo, 
			       (PLA_Obj) *alpha, (PLA_Obj) *x, (PLA_Obj) *y, 
			       (PLA_Obj) *a );

  return;
}

#if MANUFACTURE == CRAY
void PLA_TRMV_F
#else
void PLA_FOR2C( pla_trmv_f )
#endif
( FINT *uplo, FINT *trans, FINT *diag, FINT *a, FINT *x, 
		 FINT *ierror )
{
  *ierror = ( FINT ) PLA_Trmv( ( int ) *uplo, ( int ) *trans, ( int ) *diag, 
			       (PLA_Obj) *a, (PLA_Obj) *x ); 

  return;
}


#if MANUFACTURE == CRAY
void PLA_TRSV_F
#else
void PLA_FOR2C( pla_trsv_f )
#endif
( FINT *uplo, FINT *trans, FINT *diag, FINT *a, FINT *x, 
		 FINT *ierror )
{
  *ierror = ( FINT ) PLA_Trsv( ( int ) *uplo, ( int ) *trans, ( int ) *diag, 
			       (PLA_Obj) *a, (PLA_Obj) *x ); 

  return;
}

/*************************************************************
   Chapter 8
**************************************************************/

/*
   Local level-3 BLAS
*/

#if MANUFACTURE == CRAY
void PLA_LOCALGEMM_F
#else
void PLA_FOR2C( pla_localgemm_f )
#endif
( FINT *transa, FINT *transb, 
		       FINT *alpha, FINT *a, FINT *b, 
		       FINT *beta,  FINT *c, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_gemm( ( int ) *transa, ( int ) *transb, 
			(PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *b, 
		        (PLA_Obj) *beta,  (PLA_Obj) *c );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALSYMM_F
#else
void PLA_FOR2C( pla_localsymm_f )
#endif
( FINT *side, FINT *uplo, FINT *alpha, FINT *a, FINT *b, 
		       FINT *beta,  FINT *c, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_symm( ( int ) *side, ( int ) *uplo,
			(PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *b, 
		        (PLA_Obj) *beta,  (PLA_Obj) *c );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALSYRK_F
#else
void PLA_FOR2C( pla_localsyrk_f )
#endif
( FINT *uplo, FINT *trans, 
		      FINT *alpha, FINT *a, 
                      FINT *beta,  FINT *c, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_syrk( ( int ) *uplo, ( int ) *trans,
			(PLA_Obj) *alpha, (PLA_Obj) *a, 
		        (PLA_Obj) *beta,  (PLA_Obj) *c );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALSYR2K_F
#else
void PLA_FOR2C( pla_localsyr2k_f )
#endif
( FINT *uplo, FINT *trans, FINT *alpha, FINT *a, FINT *b,
			FINT *beta,  FINT *c, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_syr2k( ( int ) *uplo, ( int ) *trans,
			(PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *b, 
		        (PLA_Obj) *beta,  (PLA_Obj) *c );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALTRMM_F
#else
void PLA_FOR2C( pla_localtrmm_f )
#endif
( FINT *side, FINT *uplo, FINT *trans, FINT *diag,
		       FINT *alpha, FINT *a, FINT *b, FINT *ierror ) 
{
  *ierror = ( FINT ) PLA_Local_trmm( ( int ) *side,  ( int ) *uplo, 
				     ( int ) *trans, ( int ) *diag,
				     (PLA_Obj) *alpha, (PLA_Obj) *a, 
				     (PLA_Obj) *b );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LOCALTRSM_F
#else
void PLA_FOR2C( pla_localtrsm_f )
#endif
( FINT *side, FINT *uplo, FINT *trans, FINT *diag,
		       FINT *alpha, FINT *a, FINT *b, FINT *ierror ) 
{
  *ierror = ( FINT ) PLA_Local_trsm( ( int ) *side,  ( int ) *uplo, 
				     ( int ) *trans, ( int ) *diag,
				     (PLA_Obj) *alpha, (PLA_Obj) *a, 
				     (PLA_Obj) *b );

  return;
}


/*
   Global level-3 BLAS
*/

#if MANUFACTURE == CRAY
void PLA_GEMM_F
#else
void PLA_FOR2C( pla_gemm_f )
#endif
( FINT *transa, FINT *transb, FINT *alpha, FINT *a, FINT *b, 
		 FINT *beta,  FINT *c, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Gemm( ( int ) *transa, ( int ) *transb, 
			(PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *b, 
		        (PLA_Obj) *beta,  (PLA_Obj) *c );

  return;
}

#if MANUFACTURE == CRAY
void PLA_SYMM_F
#else
void PLA_FOR2C( pla_symm_f )
#endif
( FINT *side, FINT *uplo, FINT *alpha, FINT *a, FINT *b, 
		 FINT *beta,  FINT *c, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Symm( ( int ) *side, ( int ) *uplo,
			       (PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *b, 
			       (PLA_Obj) *beta,  (PLA_Obj) *c );

  return;
}

#if MANUFACTURE == CRAY
void PLA_SYRK_F
#else
void PLA_FOR2C( pla_syrk_f )
#endif
( FINT *uplo, FINT *trans, 
		 FINT *alpha, FINT *a, 
		 FINT *beta,  FINT *c, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Syrk( ( int ) *uplo, ( int ) *trans,
			       (PLA_Obj) *alpha, (PLA_Obj) *a, 
			       (PLA_Obj) *beta,  (PLA_Obj) *c );

  return;
}

#if MANUFACTURE == CRAY
void PLA_SYR2K_F
#else
void PLA_FOR2C( pla_syr2k_f )
#endif
( FINT *uplo, FINT *trans, FINT *alpha, FINT *a, FINT *b,
		  FINT *beta,  FINT *c, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Syr2k( ( int ) *uplo, ( int ) *trans,
				(PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *b, 
				(PLA_Obj) *beta,  (PLA_Obj) *c );

  return;
}

#if MANUFACTURE == CRAY
void PLA_TRMM_F
#else
void PLA_FOR2C( pla_trmm_f )
#endif
( FINT *side, FINT *uplo, FINT *trans, FINT *diag,
		 FINT *alpha, FINT *a, FINT *b, FINT *ierror ) 
{
  *ierror = ( FINT ) PLA_Trmm( ( int ) *side,  ( int ) *uplo, 
			       ( int ) *trans, ( int ) *diag,
			       (PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *b );

  return;
}

#if MANUFACTURE == CRAY
void PLA_TRSM_F
#else
void PLA_FOR2C( pla_trsm_f )
#endif
( FINT *side, FINT *uplo, FINT *trans, FINT *diag,
		 FINT *alpha, FINT *a, FINT *b, FINT *ierror ) 
{
  *ierror = ( FINT ) PLA_Trsm( ( int ) *side,  ( int ) *uplo,
			       ( int ) *trans, ( int ) *diag,
			       (PLA_Obj) *alpha, (PLA_Obj) *a, (PLA_Obj) *b );

  return;
}

#if MANUFACTURE == CRAY
void PLA_ENVIRONNBALG_F
#else
void PLA_FOR2C( pla_environnbalg_f )
#endif
( FINT *operation, FINT *template, FINT *nb_alg, 
			   FINT *ierror )
{
  int
    nb_alg_C;

  *ierror = ( FINT ) PLA_Environ_nb_alg( ( int) *operation, (PLA_Template) *template, &nb_alg_C );

  *nb_alg = ( FINT ) nb_alg_C;

  return;
}


/*************************************************************
   Miscellaneous
*************************************************************/

/*
#if MANUFACTURE == CRAY
void PLA_COPYADD_F
#else
void PLA_FOR2C( pla_copyadd_f )
#endif
( FINT *obj_from, FINT *obj_to, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Copy_add( (PLA_Obj) *obj_from, (PLA_Obj) *obj_to );

  return;
}
*/

/*
#if MANUFACTURE == CRAY
void PLA_DOTADD_F
#else
void PLA_FOR2C( pla_dotadd_f )
#endif
( FINT *x, FINT *y, FINT *alpha , FINT *ierror )
{
  *ierror = ( FINT ) PLA_Dot_add( (PLA_Obj) *x, (PLA_Obj) *y, (PLA_Obj) *alpha );

  return;
}
*/

/*
#if MANUFACTURE == CRAY
void PLA_DOTX_F
#else
void PLA_FOR2C( pla_dotx_f )
#endif
( FINT *alpha, FINT *x, FINT *y, FINT *beta, FINT *result,
                  FINT *ierror )
{
  *ierror = ( FINT ) PLA_Dot_x( (PLA_Obj) *alpha, (PLA_Obj) *x, (PLA_Obj) *y, 
				(PLA_Obj) *beta, (PLA_Obj) *result );

  return;
}
*/

#if MANUFACTURE == CRAY
void PLA_CREATECONSTANTSCONFTO_F
#else
void PLA_FOR2C( pla_createconstantsconfto_f )
#endif
( FINT *a, 
				     FINT *minus_one, 
				     FINT *zero,
				     FINT *one, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Create_constants_conf_to( (PLA_Obj) *a,
						   (PLA_Obj *) minus_one,
						   (PLA_Obj *) zero,
						   (PLA_Obj *) one );

  return;
}


#if MANUFACTURE == CRAY
void PLA_LOCALEQUALZERO_F
#else
void PLA_FOR2C( pla_localequalzero_f )
#endif
( FINT *alpha, FINT *result, FINT *ierror )
{
  *ierror = ( FINT ) PLA_SUCCESS;
  *result = ( FINT ) PLA_Local_equal_zero( (PLA_Obj) *alpha );

  return;
}


#if MANUFACTURE == CRAY
void PLA_LOCALINVERTSIGN_F
#else
void PLA_FOR2C( pla_localinvertsign_f )
#endif
( FINT *alpha, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_invert_sign( (PLA_Obj) *alpha );

  return; 
}


#if MANUFACTURE == CRAY
void PLA_LOCALSIGN_F
#else
void PLA_FOR2C( pla_localsign_f )
#endif
( FINT *alpha, FINT *sign, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_sign( (PLA_Obj) *alpha, (PLA_Obj) *sign );

  return; 
}


/*
#if MANUFACTURE == CRAY
void PLA_LOCALSQRT_F
#else
void PLA_FOR2C( pla_localsqrt_f )
#endif
( FINT *alpha, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_sqrt( (PLA_Obj) *alpha ); 

  return; 
}
*/


#if MANUFACTURE == CRAY
void PLA_OBJSETTOIDENTITY_F
#else
void PLA_FOR2C( pla_objsettoidentity_f )
#endif
( FINT *a, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Obj_set_to_identity( (PLA_Obj) *a );

  return; 
}


#if MANUFACTURE == CRAY
void PLA_OBJSPLITSIZETONEXTPROC_F
#else
void PLA_FOR2C( pla_objsplitsizetonextproc_f )
#endif
( FINT *obj, FINT *side, FINT *size, 
					FINT *owner, FINT *ierror )
{
  int
    size_C, owner_C;

  *ierror = ( FINT ) PLA_Obj_split_size_to_next_proc 
				   ( (PLA_Obj) *obj, ( int ) *side, 
				     &size_C, &owner_C );

  *size = ( FINT ) size_C;
  *owner = ( FINT ) owner_C;

  return;
}


/*
#if MANUFACTURE == CRAY
void PLA_SCALX_F
#else
void PLA_FOR2C( pla_scalx_f )
#endif
( FINT *shape, FINT *alpha, FINT *x, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Scal_x( ( int ) *shape, 
				 (PLA_Obj) *alpha, (PLA_Obj) *x );

  return;
}
*/

#if MANUFACTURE == CRAY
void PLA_SHOW_F
#else
void PLA_FOR2C( pla_show_f )
#endif
( FINT *obj, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Show( (PLA_Obj) *obj );

  return;
}


#if MANUFACTURE == CRAY
void PLA_SYRKPERFORMLOCALPART_F
#else
void PLA_FOR2C( pla_syrkperformlocalpart_f )
#endif
( FINT *uplo, 
		      FINT *alpha, FINT *a, 
                      FINT *beta,  FINT *c, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Syrk_perform_local_part( ( int ) *uplo,  
			(PLA_Obj) *alpha, (PLA_Obj) *a, 
		        (PLA_Obj) *beta,  (PLA_Obj) *c );

  return;
}

#if MANUFACTURE == CRAY
void PLA_LU_F
#else
void PLA_FOR2C( pla_lu_f )
#endif
( FINT *A, FINT *pivots, FINT *ierror )
{
  *ierror = ( FINT ) PLA_LU( (PLA_Obj) *A, (PLA_Obj) *pivots );

  return;
}

#if MANUFACTURE == CRAY
void PLA_APPLYPIVOTSTOROWS_F
#else
void PLA_FOR2C( pla_applypivotstorows_f )
#endif
( FINT *A, FINT *pivots, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Apply_pivots_to_rows( (PLA_Obj) *A, (PLA_Obj) *pivots );

  return;
}

#if MANUFACTURE == CRAY
void PLA_PIVOTFROMRIGHTREV_F
#else
void PLA_FOR2C( pla_pivotfromrightrev_f )
#endif
( FINT *A, FINT *pivots, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Apply_pivots_to_columns_in_reverse( (PLA_Obj) *A, (PLA_Obj) *pivots );

  return;
}


#if MANUFACTURE == CRAY
void PLA_CHOL_F
#else
void PLA_FOR2C( pla_chol_f )
#endif
( FINT *uplo, FINT *A, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Chol( (int) *uplo, (PLA_Obj) *A );

  return;
}


#if MANUFACTURE == CRAY
void PLA_LOCALCHOL_F
#else
void PLA_FOR2C( pla_localchol_f )
#endif
( FINT *uplo, FINT *A, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Local_chol( (int) *uplo, (PLA_Obj) *A );

  return;
}

#if MANUFACTURE == CRAY
void PLA_APPLYWYTRANSFORM_F
#else
void PLA_FOR2C( pla_applywytransform_f )
#endif
( FINT *uplo, FINT *trans, FINT *W, FINT *Y, FINT *A, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Apply_W_Y_transform( 
	   (int) *uplo, (int) *trans, (PLA_Obj) *W, (PLA_Obj) *Y, 
	   (PLA_Obj) *A );

  return;
}

#if MANUFACTURE == CRAY
void PLA_COMPUTEHOUSEV_F
#else
void PLA_FOR2C( pla_computehousev_f )
#endif
( FINT *x, FINT *beta, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Compute_House_v( (PLA_Obj) *x, (PLA_Obj) *beta ); 

  return;
}

#if MANUFACTURE == CRAY
void PLA_COMPUTEWY_F
#else
void PLA_FOR2C( pla_computewy_f )
#endif
( FINT *A_mv, FINT *s, FINT *W_mv, FINT *Y_mv, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Compute_WY( (PLA_Obj) *A_mv, (PLA_Obj) *s,
				      (PLA_Obj) *W_mv, (PLA_Obj) *Y_mv ); 

  return;
}

#if MANUFACTURE == CRAY
void PLA_FORMQ_F
#else
void PLA_FOR2C( pla_formq_f )
#endif
( FINT *trans, FINT *A, FINT *s, FINT *Q, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Form_Q( (int) *trans, (PLA_Obj) *A, (PLA_Obj) *s,
				      (PLA_Obj) *Q );

  return;
}

#if MANUFACTURE == CRAY
void PLA_QR_F
#else
void PLA_FOR2C( pla_qr_f )
#endif
( FINT *A, FINT *s, FINT *ierror )
{
  *ierror = ( FINT ) PLA_QR( (PLA_Obj) *A, (PLA_Obj) *s );

  return;
}


#if MANUFACTURE == CRAY
void PLA_QSOLVE_F
#else
void PLA_FOR2C( pla_q_solve_f )
#endif
( FINT *side, FINT *trans, FINT *A, FINT *s, FINT *B, FINT *ierror )
{
  *ierror = ( FINT ) PLA_Q_solve( (int) *side, (int) *trans, (PLA_Obj) *A, (PLA_Obj) *s,
				      (PLA_Obj) *B );

  return;
}



#if MANUFACTURE == CRAY
void PLA_SVD_F
#else
void PLA_FOR2C( pla_svd_f )
#endif
( FINT *A, FINT *U, FINT *D, FINT *V, FINT *ierror )
{
  *ierror = ( FINT ) PLA_SVD( (PLA_Obj) *A, (PLA_Obj) *U,  (PLA_Obj) *D, 
				      (PLA_Obj) *V );

  return;
}








