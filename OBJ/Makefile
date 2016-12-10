include ../Make.include

OBJ = 	PLA_Obj.o \
	PLA_Obj_enter_exit.o \
	PLA_Obj_management.o \
	PLA_Memory.o \
	PLA_Diagnostics.o \
	PLA_Utility.o \
	PLA_Temp_extra_routines.o \
	PLA_Temp_extra_routines_enter_exit.o \
	PLA_Obj_extra_routines.o \
	PLA_Obj_set_to_identity.o \
	PLA_Obj_ooc.o \
	PLA_Obj_ooc_enter_exit.o \
	PLA_Timings.o \
	PLA_Check_R12.o 

OBJ: $(OBJ)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch OBJ

clean:
	rm -f $(OBJ) core *.o *~

