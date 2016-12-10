include ../Make.include

POOCLAPACK: 
	cd OOC_BLAS3; make OOC_BLAS3
	cd OOC_Cholesky; make OOC_Cholesky
	rm -f POOCLAPACK

#	$(AR) vr $(PLAPACKLIB) 
#	$(RANLIB) $(PLAPACKLIB)
	touch POOCLAPACK

clean:
	cd OOC_BLAS3; make clean
	cd OOC_Cholesky; make clean
	cd TEST; make clean
	rm -f POOCLAPACK

