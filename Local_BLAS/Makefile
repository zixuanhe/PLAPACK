include ../Make.include

local_blas = 	\
        PLA_Local_BLAS1.o \
	PLA_Local_BLAS2.o \
	PLA_Local_BLAS1_enter_exit.o \
	PLA_Local_BLAS2_enter_exit.o \
	PLA_Local_BLAS3_enter_exit.o \
	PLA_Local_iamax.o \
	PLA_Local_trmv.o \
	PLA_Local_syr.o \
	PLA_Local_syr2.o \
        PLA_Local_gemm.o \
        PLA_Local_symm.o \
	PLA_Local_syrk.o \
	PLA_Local_herk.o \
	PLA_Local_syr2k.o \
	PLA_Local_trsm.o \
	PLA_Local_trmm.o \
	PLA_Local_inquiry.o \
	PLA_Local_misc.o \
	PLA_Local_absolute_max.o \
	PLA_Conjugate.o \
	PLA_OMP_gemm.o

#	PLA_Local_trsv.o 
#        PLA_Local_gemv.o 

Local_BLAS: $(local_blas)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch Local_BLAS

clean:
	rm -f $(local_blas) core *.o *~

