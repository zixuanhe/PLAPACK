include ../Make.include

blas = 	PLA_Gemm.o \
	PLA_Gemm_enter_exit.o \
        PLA_Gemm_x.o \
	PLA_Gemm_x_enter_exit.o \
	PLA_Gemm_panpan.o \
	PLA_Gemm_matpan.o \
	PLA_Gemm_panmat.o \
	PLA_Gemm_A.o \
	PLA_Gemm_B.o \
	PLA_Gemm_C.o \
	PLA_Symm.o \
	PLA_Symm_enter_exit.o \
	PLA_Trmm.o \
	PLA_Trmm_enter_exit.o \
	PLA_Trmm_left_lower.o \
	PLA_Trmm_left_upper.o \
	PLA_Trmm_right_lower.o \
	PLA_Trmm_right_upper.o \
	PLA_Trsm_right_lower.o \
	PLA_Trsm.o \
	PLA_Trsm_enter_exit.o \
	PLA_Trsm_left_lower.o \
	PLA_Trsm_left_upper.o \
	PLA_Trsm_right_lower.o \
	PLA_Trsm_right_upper.o \
	PLA_Syrk.o \
	PLA_Syrk_enter_exit.o \
	PLA_Syrk_panpan.o \
	PLA_Syrk_perform_local_part.o \
	PLA_Herk_perform_local_part.o \
	PLA_Syr2k.o 

BLAS3: $(blas)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch BLAS3

clean:
	rm -f $(blas) core *.o *~

