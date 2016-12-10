include ../../Make.include

OBJS = 	PLA_Local_bidiag_svd.o \
	PLA_SVD.o  \
	PLA_SVD_enter_exit.o  \
	dbdsqr_x.o

SVD:    $(OBJS)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch SVD
clean:
	rm -f $(OBJS) SVD core *.o *~

