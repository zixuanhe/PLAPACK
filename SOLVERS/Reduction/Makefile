include ../../Make.include

OBJS = 	PLA_Apply_sym_House.o \
	PLA_Tri_red.o \
	PLA_Tri_red_enter_exit.o \
	PLA_Create_sym_eigenproblem.o \
	PLA_Bi_red_simple.o \
	PLA_Bi_red_enter_exit.o \
	PLA_Create_svd_problem.o


Red: $(OBJS)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch Red
clean:
	rm -f $(OBJS) Red core *.o *~

