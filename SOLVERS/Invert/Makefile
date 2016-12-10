include ../../Make.include

OBJS = 	PLA_General_invert.o \
	PLA_General_invert_factors.o \
	PLA_General_invert_inv.o \
	PLA_General_invert_stable.o \
	PLA_General_invert_enter_exit.o \
	PLA_Triangular_invert.o \
	PLA_Triangular_invert_lower.o \
	PLA_Triangular_invert_enter_exit.o \
	PLA_Triangular_At_A.o \
	PLA_Triangular_At_A_enter_exit.o \
	PLA_Triangular_Lt_L.o \
	PLA_Symmetric_invert.o \
	PLA_Symmetric_invert_lower.o \
	PLA_Symmetric_invert_enter_exit.o \
	PLA_Local_tri_inv.o

Inv: $(OBJS)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch Inv
clean:
	rm -f $(OBJS) Inv core *.o *~

