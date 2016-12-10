PLAPACK_ROOT = ../..
include $(PLAPACK_ROOT)/Make.include

TMP = /tmp

UTIL   = create_problem.o \
	create_problem_ooc.o \
	create_problem_ooc_spd.o \
	create_problem_ooc_by_panels_spd.o

TEST_OOC_SYRK = test_ooc_syrk.o

TEST_OOC_SYRK_ASYNC = test_ooc_syrk_async.o

TEST_OOC_LU = test_ooc_lu.o 

TEST_OOC_CHOL = test_ooc_chol.o 

TEST_OOC_CHOL_ASYNC = test_ooc_chol_async.o

TEST_OOC_CHOL_TILES = test_ooc_chol_by_panels.o 

TEST_OOC_CHOL_ASYNC_TILES = test_ooc_chol_async_by_panels.o 

wcr:
	make test_chol
	make test_ooc_chol
	make test_ooc_chol_async
	make test_ooc_chol_by_panels
	make test_ooc_chol_async_by_panels

test_ooc_lu: $(TEST_OOC_LU) $(UTIL)
	rm -f TEST_OOC_LU
	$(LINKER) -o $(TMP)/test_ooc_lu.x $(LFLAGS) $(TEST_OOC_LU) $(UTIL) \
        $(PLAPACKLIB) $(PLAPACKLIB) $(LIB)
	chmod a+rx $(TMP)/test_ooc_lu.x
	touch TEST_OOC_LU

test_ooc_syrk: $(TEST_OOC_SYRK) $(UTIL)
	rm -f TEST_OOC_SYRK
	$(LINKER) -o $(TMP)/test_ooc_syrk.x $(LFLAGS) $(TEST_OOC_SYRK) $(UTIL)\
        $(PLAPACKLIB) $(PLAPACKLIB) $(LIB)
	touch TEST_OOC_SYRK
	chmod a+rx $(TMP)/test_ooc_syrk.x

test_ooc_syrk_async: $(TEST_OOC_SYRK_ASYNC) $(UTIL)
	rm -f TEST_OOC_SYRK_ASYNC
	$(LINKER) -o $(TMP)/test_ooc_syrk_async.x $(LFLAGS) $(TEST_OOC_SYRK_ASYNC) $(UTIL) \
        $(PLAPACKLIB) $(PLAPACKLIB) $(LIB)
	touch TEST_OOC_SYRK_ASYNC
	chmod a+rx $(TMP)/test_ooc_syrk_async.x

test_ooc_chol: $(TEST_OOC_CHOL) $(UTIL)
	rm -f TEST_OOC_CHOL
	$(LINKER) -o $(TMP)/test_ooc_chol.x $(LFLAGS) $(TEST_OOC_CHOL) $(UTIL)\
        $(PLAPACKLIB) $(PLAPACKLIB) $(LIB)
	chmod a+rx $(TMP)/test_ooc_chol.x
	touch TEST_OOC_CHOL

test_ooc_chol_async: $(TEST_OOC_CHOL_ASYNC) $(UTIL)
	rm -f TEST_OOC_CHOL_ASYNC
	$(LINKER) -o $(TMP)/test_ooc_chol_async.x $(LFLAGS) $(TEST_OOC_CHOL_ASYNC) $(UTIL) \
        $(PLAPACKLIB) $(PLAPACKLIB) $(LIB)
	chmod a+rx $(TMP)/test_ooc_chol_async.x
	touch TEST_OOC_CHOL_ASYNC

test_ooc_chol_by_panels: $(TEST_OOC_CHOL_TILES) $(UTIL)
	rm -f TEST_OOC_CHOL_TILES
	$(LINKER) -o $(TMP)/test_ooc_chol_by_panels.x $(LFLAGS) $(TEST_OOC_CHOL_TILES) $(UTIL) \
        $(PLAPACKLIB) $(PLAPACKLIB) $(LIB)
	chmod a+rx $(TMP)/test_ooc_chol_by_panels.x
	touch TEST_OOC_CHOL_TILES

test_ooc_chol_async_by_panels: $(TEST_OOC_CHOL_ASYNC_TILES) $(UTIL)
	rm -f TEST_OOC_CHOL_ASYNC_TILES
	$(LINKER) -o $(TMP)/test_ooc_chol_async_by_panels.x $(LFLAGS) $(TEST_OOC_CHOL_ASYNC_TILES) $(UTIL) \
        $(PLAPACKLIB) $(PLAPACKLIB) $(LIB)
	chmod a+rx $(TMP)/test_ooc_chol_async_by_panels.x
	touch TEST_OOC_CHOL_ASYNC_TILES

clean :
	rm -f TEST_* *.o core *~ file* *.x

