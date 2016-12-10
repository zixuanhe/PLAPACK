/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

struct PLA_Template_struct{
  MPI_Comm
    comm_all,
    comm_row,
    comm_col;

  int
    nb, 
    zero_or_one,

    me,
    myrow,
    mycol,
    
    nprows,
    npcols,
    
    times_used;
}; 

typedef struct PLA_Template_struct *PLA_Template;
