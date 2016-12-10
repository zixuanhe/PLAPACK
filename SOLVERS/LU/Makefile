include ../../Make.include

OBJS = 	PLA_LU.o \
	PLA_LU_enter_exit.o \
	PLA_LU_mv.o \
	PLA_LU_mv_enter_exit.o \
	PLA_Apply_pivots_new.o\
	PLA_seq_lu.o \
	PLA_Local_lu.o \
	PLA_General.o

LU: $(OBJS)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch LU
clean:
	rm -f $(OBJS) LU core *.o *~

