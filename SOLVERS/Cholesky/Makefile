include ../../Make.include

OBJS = 	PLA_Chol.o \
	PLA_Chol_simple.o \
	PLA_Local_chol.o \
	PLA_Chol_enter_exit.o \
	PLA_Pos_def.o \
	PLA_Pos_def_enter_exit.o


Cholesky: $(OBJS)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch Cholesky
clean:
	rm -f $(OBJS) Cholesky core *.o *~

