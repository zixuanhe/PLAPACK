PLAPACK_ROOT = ..

include $(PLAPACK_ROOT)/Make.include

UTIL   = create_problem.o create_problem_mv.o

UTIL_OOC = create_problem_ooc.o $(UTIL)

UTIL_OLD = create_problem_old.o PLA_Gemm_panpan.o

TEST_LU = test_lu.o 

TEST_ELEM_MULT = test_elem_mult.o 

TEST_INV = test_inv.o 

TEST_TRI_INV = test_tri_inv.o 

TEST_SYM_INV = test_sym_inv.o 

TEST_SPECTRAL_DECOMP = test_spectral_decomp.o 

TEST_SVD = test_svd.o 

TEST_TRI_RED = test_tri_red.o 

TEST_SHIFT = test_shift.o

TEST_MATRIX_ONE_NORM = test_matrix_one_norm.o

TEST_MATRIX_INF_NORM = test_matrix_inf_norm.o

TEST_IAMAX = test_iamax.o

TEST_SYRK = test_syrk.o 

TEST_OOC_SYRK = test_ooc_syrk.o 

TEST_GEMM = test_gemm.o 

TEST_TRMM = test_trmm.o 

TEST_TRSM = test_trsm.o 

TEST_CHOL = test_chol.o 

TEST_GEMV = test_gemv.o 

TEST_TRMV = test_trmv.o 

TEST_DIAG = test_diag.o 

test_svd: $(TEST_SVD) $(UTIL)
	rm -f TEST_SVD
	$(LINKER) -o /tmp/test_svd.x $(LFLAGS) $(TEST_SVD) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_SVD

test_elem_mult: $(TEST_ELEM_MULT) $(UTIL)
	rm -f TEST_ELEM_MULT
	$(LINKER) -o /tmp/test_elem_mult.x $(LFLAGS) $(TEST_ELEM_MULT) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_ELEM_MULT

test_diag: $(TEST_DIAG) $(UTIL)
	rm -f TEST_DIAG
	$(LINKER) -o /tmp/test_diag.x $(LFLAGS) $(TEST_DIAG) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_DIAG


test_sym_inv: $(TEST_SYM_INV) $(UTIL)
	rm -f TEST_SYM_INV
	$(LINKER) -o /tmp/test_sym_inv.x $(LFLAGS) $(TEST_SYM_INV) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_SYM_INV

test_tri_inv: $(TEST_TRI_INV) $(UTIL)
	rm -f TEST_TRI_INV
	$(LINKER) -o /tmp/test_tri_inv.x $(LFLAGS) $(TEST_TRI_INV) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_TRI_INV

test_inv: $(TEST_INV) $(UTIL)
	rm -f TEST_INV
	$(LINKER) -o /tmp/test_inv.x $(LFLAGS) $(TEST_INV) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_INV

test_spectral_decomp: $(TEST_SPECTRAL_DECOMP) $(UTIL)
	rm -f TEST_SPECTRAL_DECOM
	$(LINKER) -o /tmp/test_spectral_decomp.x $(LFLAGS) $(TEST_SPECTRAL_DECOMP) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_SPECTRAL_DECOMP

test_tri_red: $(TEST_TRI_RED) $(UTIL)
	rm -f TEST_TRI_RED
	$(LINKER) -o /tmp/test_tri_red.x $(LFLAGS) $(TEST_TRI_RED) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_TRI_RED

test_matrix_one_norm: $(TEST_MATRIX_ONE_NORM) $(UTIL)
	rm -f TEST_MATRIX_ONE_NORM
	$(LINKER) -o /tmp/test_matrix_one_norm.x $(LFLAGS) $(TEST_MATRIX_ONE_NORM) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_MATRIX_ONE_NORM

test_matrix_inf_norm: $(TEST_MATRIX_INF_NORM) $(UTIL)
	rm -f TEST_MATRIX_INF_NORM
	$(LINKER) -o /tmp/test_matrix_inf_norm.x $(LFLAGS) $(TEST_MATRIX_INF_NORM) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_MATRIX_INF_NORM

test_trsm: $(TEST_TRSM) $(UTIL)
	rm -f TEST_TRSM
	$(LINKER) -o /tmp/test_trsm.x $(LFLAGS) $(TEST_TRSM) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_TRSM

test_lu: $(TEST_LU) $(UTIL)
	rm -f TEST_LU
	$(LINKER) -o /tmp/test_lu.x $(LFLAGS) $(TEST_LU) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_LU

test_shift: $(TEST_SHIFT) $(UTIL)
	rm -f TEST_SHIFT
	$(LINKER) -o /tmp/test_shift.x $(LFLAGS) $(TEST_SHIFT) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_SHIFT

test_iamax: $(TEST_IAMAX) $(UTIL)
	rm -f TEST_IAMAX
	$(LINKER) -o /tmp/test_iamax.x $(LFLAGS) $(TEST_IAMAX) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_IAMAX

test_ooc_syrk: $(TEST_OOC_SYRK) $(UTIL_OOC)
	rm -f TEST_OOC_SYRK
	$(LINKER) -o /tmp/test_ooc_syrk.x $(LFLAGS) $(TEST_OOC_SYRK) $(UTIL_OOC) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_OOC_SYRK

test_syrk: $(TEST_SYRK) $(UTIL)
	rm -f TEST_SYRK
	$(LINKER) -o /tmp/test_syrk.x $(LFLAGS) $(TEST_SYRK) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_SYRK

test_ooc: $(TEST_OOC) $(UTIL) 
	rm -f TEST_OOC
	$(LINKER) -o /tmp/test_ooc.x $(LFLAGS) $(TEST_OOC) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_OOC

test_gemm: $(TEST_GEMM) $(UTIL)
	rm -f TEST_GEMM
	$(LINKER) -o /tmp/test_gemm.x $(LFLAGS) $(TEST_GEMM) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_GEMM

test_trmm: $(TEST_TRMM) $(UTIL)
	rm -f TEST_TRMM
	$(LINKER) -o /tmp/test_trmm.x $(LFLAGS) $(TEST_TRMM) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_TRMM

test_chol: $(TEST_CHOL) $(UTIL)
	rm -f TEST_CHOL
	$(LINKER) -o /tmp/test_chol.x $(LFLAGS) $(TEST_CHOL) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_CHOL

test_gemv: $(TEST_GEMV) $(UTIL)
	rm -f TEST_GEMV
	$(LINKER) -o /tmp/test_gemv.x $(LFLAGS) $(TEST_GEMV) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_GEMV

test_trmv: $(TEST_TRMV) $(UTIL)
	rm -f TEST_TRMV
	$(LINKER) -o /tmp/test_trmv.x $(LFLAGS) $(TEST_TRMV) $(UTIL) $(OBJS) \
        $(PLAPACKLIB) $(MMPILIB) $(lapackLIB) $(LIB)
	touch TEST_TRMV

clean :
	rm -f TEST_* *.o MMPI* core *~ file* *.x output_test_gemm *log

