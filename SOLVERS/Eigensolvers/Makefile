include ../../Make.include

OBJS = 	PLA_Local_sym_tridiag_eig.o \
	PLA_Spectral_decomp.o \
	dsteqr_x.o \
	ssteqr_x.o \
	ssyev.o 

Eig:    $(OBJS)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch Eig
clean:
	rm -f $(OBJS) Eig core *.o *~

