include ../../Make.include

OBJS = \
                PLA_QR.o \
                PLA_QR_enter_exit.o \
		PLA_QR_right.o \
		PLA_QR_mv.o \
		PLA_Compute_House_v.o \
		PLA_Apply_House.o \
		PLA_Compute_W_Y.o \
		PLA_Apply_W_Y.o \
	        PLA_Form_Q.o  \
	        PLA_Form_Q_x.o  \
		PLA_Q_solve.o \
	 	PLA_Q_solve_mv.o \
	 	PLA_Q_solve_matrix.o \
	        PLA_QR_update.o \
		PLA_Least_sqr_solve.o


QR: $(OBJS)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch QR

clean :
	rm -f $(OBJS) QR core *.o *~ 

