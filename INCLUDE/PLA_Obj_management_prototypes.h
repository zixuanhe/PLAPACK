/*
   PLAPACK Release 3.0
   
   Copyright (C) 2000
   Robert A. van de Geijn and The University of Texas at Austin

   For GNU-license details see the file GNU_license in the 
   PLAPACK root directory.
*/

int PLA_Obj_attach_buffer ( PLA_Obj, void *, int, int );

/*----------------------------------------------------------------------------

Purpose : Attach a given data buffer to the object.

IN     obj               object to which buffer is to be attached
IN     buffer            buffer to be used for data
IN     ldim              leading dimension of buffer
IN     user_buffer       indicates whether PLAPACK provided buffer,
                         or user provides buffer.

----------------------------------------------------------------------------*/

PLA_Obj_base pla_get_base_obj ( );

PLA_Obj pla_get_view ( );

int pla_free_view( PLA_Obj );

int pla_detach_base_obj( PLA_Obj_base );

