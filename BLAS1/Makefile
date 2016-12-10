include ../Make.include

blas = 	PLA_Axpy.o \
	PLA_Axpy_enter_exit.o \
	PLA_Swap.o \
        PLA_Iamax.o \
        PLA_Iamax_enter_exit.o \
	PLA_Inv_scal.o \
	PLA_Inv_scal_enter_exit.o \
	PLA_Scal.o \
	PLA_Scal_enter_exit.o \
        PLA_Asum.o \
        PLA_Asum_enter_exit.o \
        PLA_Nrm2.o \
        PLA_Nrm2_enter_exit.o \
        PLA_Dot.o \
        PLA_Dot_enter_exit.o

BLAS1: $(blas)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch BLAS1

clean:
	rm -f $(blas) *.o core *~

