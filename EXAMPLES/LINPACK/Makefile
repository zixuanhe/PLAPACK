include ../../Make.include

OBJS = create_problem.o driver.o process_answer.o \
       PLA_LU_linpack.o

LINPACK : $(OBJS) $(PLAPACKLIB) 
	rm -f LU_driver*
	$(LINKER) -o /tmp/linpack.x $(LFLAGS) $(OBJS) $(PLAPACKLIB) $(LIB)
	touch LU_driver

clean :
	rm -f LU_driver* *.o core *~

