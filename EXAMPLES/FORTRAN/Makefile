include ../../Make.include

OBJS_LU   = create_problem_f.o driver.o process_answer_f.o 

LU_driver : $(OBJS_LU) $(PLAPACKLIB)
	rm -f LU_driver*
	$(FC) -o /tmp/LU_driver.x $(FLFLAGS) $(OBJS_LU) $(PLAPACKLIB) $(LIB)
	touch LU_driver

clean :
	rm -f LU_driver* *.o core *~

