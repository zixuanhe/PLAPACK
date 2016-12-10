include ../Make.include

COPY =	PLA_Copy.o \
	PLA_Copy_enter_exit.o \
	PLA_Copy_util.o \
	PLA_Pack_unpack.o \
	PLA_Compute_subbuffer_info.o \
	PLA_Copy_x.o \
	PLA_OOC_Copy.o \
	PLA_Copy_sym_tridiag_to_msc.o \
	PLA_Copy_bidiag_to_msc.o \
	PLA_Copy_diag_to_mv.o \
	PLA_MPI.o

COPY: $(COPY)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch COPY

clean:
	rm -f $(COPY) *~ core

