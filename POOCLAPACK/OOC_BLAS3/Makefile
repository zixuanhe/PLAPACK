include ../../Make.include

ooc_blas3 = \
	PLA_OOC_Syrk_enter_exit.o \
	PLA_OOC_Syrk.o \
	PLA_OOC_Syrk_async.o \
	PLA_OOC_Gemm.o \
	PLA_OOC_Gemm_async.o 

OOC_BLAS3: $(ooc_blas3)
	$(AR) vr ../../libPLAPACK.a   $? 
	$(RANLIB) ../../libPLAPACK.a 
	touch OOC_BLAS3

clean:
	rm -f $(ooc_blas3) *~ OOC_BLAS3 core

