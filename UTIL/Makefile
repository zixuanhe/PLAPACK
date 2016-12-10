include ../Make.include

UTIL = 	PLA_Obj_split_size_to_next_proc.o \
	PLA_Global_show.o \
	PLA_Local_show.o \
	PLA_Show.o \
	PLA_Shift.o \
	PLA_Max.o \
	PLA_Set_triang_to_zero.o \
	create_random_data.o \
	pla_array_copy.o \
	PLA_Set_diagonal.o \
	PLA_Symmetrize.o

UTIL: $(UTIL)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch UTIL

clean:
	rm -f $(UTIL) *~ core

