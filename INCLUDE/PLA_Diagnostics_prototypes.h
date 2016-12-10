/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/


/***************************************************************************/

int PLA_Routine_stack_initialize( );
     
/*----------------------------------------------------------------------------

Purpose : Initialize stack that keeps track of routines called

----------------------------------------------------------------------------*/

int PLA_Routine_stack_push( char * );
     
/*----------------------------------------------------------------------------

Purpose : Push name of current routine onto routine stack

IN     routine_name      name of routine being pushed

----------------------------------------------------------------------------*/

int PLA_Routine_stack_pop( char * );
     
/*----------------------------------------------------------------------------

Purpose : Pop name of current routine from routine stack

OUT    routine_name      name of routine being popped

----------------------------------------------------------------------------*/

int PLA_Abort( char *, int, char * );
     
/*----------------------------------------------------------------------------

Purpose : Abort PLAPACK code

IN     message         nerror message

----------------------------------------------------------------------------*/

int PLA_Warning( char * );
     
/*----------------------------------------------------------------------------

Purpose : Report PLAPACK warning

IN     message         error message

----------------------------------------------------------------------------*/

int PLA_Valid_template( PLA_Template );
     
/*----------------------------------------------------------------------------

Purpose : Check if valid template

IN     templ        template to be checked

----------------------------------------------------------------------------*/

int PLA_Valid_objtype( int );
     
/*----------------------------------------------------------------------------

Purpose : Check if valid object type

IN     objtype      template to be checked

----------------------------------------------------------------------------*/
int PLA_Valid_object( PLA_Obj );
     
/*----------------------------------------------------------------------------

Purpose : Check if valid object

IN     templ        object to be checked

----------------------------------------------------------------------------*/


