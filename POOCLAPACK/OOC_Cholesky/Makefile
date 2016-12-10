include ../../Make.include

ooc_chol = \
	PLA_OOC_Chol_enter_exit.o \
	PLA_OOC_Chol.o \
	PLA_OOC_Chol_async.o \
	PLA_OOC_Chol_by_panels.o \
	PLA_OOC_Chol_async_by_panels.o 


OOC_Cholesky: $(ooc_chol)
	$(AR) vr  ../../libPLAPACK.a   $? 
	$(RANLIB) ../../libPLAPACK.a 
	touch OOC_Cholesky

clean:
	rm -f $(ooc_chol) core *~ 

