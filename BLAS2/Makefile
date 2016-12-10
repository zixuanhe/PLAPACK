include ../Make.include

blas = 	PLA_Gemv.o \
	PLA_Gemv_enter_exit.o \
	PLA_Ger.o \
	PLA_Ger_enter_exit.o \
	PLA_Syr.o \
	PLA_Syr2.o \
	PLA_Syr2_enter_exit.o \
	PLA_Syr2_perform_local_part.o \
	PLA_Symv.o \
	PLA_Symv_enter_exit.o \
	PLA_Hemv.o \
	PLA_Hemv_enter_exit.o \
	PLA_Trmv_perform_local_part.o \
	PLA_Trmv.o \
	PLA_Trmv_enter_exit.o \
	PLA_Trsv.o \
	PLA_Trsv_enter_exit.o \
	PLA_Trsv_ln.o \
	PLA_Trsv_lt.o \
	PLA_Trsv_lc.o \
	PLA_Trsv_un.o \
	PLA_Matrix_one_norm.o \
	PLA_Matrix_one_norm_enter_exit.o \
	PLA_Matrix_infinity_norm.o \
	PLA_Matrix_infinity_norm_enter_exit.o \
	PLA_Multiply_by_diagonal.o \
	PLA_Elementwise_mult.o

BLAS2: $(blas)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch BLAS2

clean:
	rm -f $(blas) core *.o *~



