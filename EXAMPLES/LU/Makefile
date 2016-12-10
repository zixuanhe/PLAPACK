include ../../Make.include

OBJS = create_problem.o driver.o process_answer.o

#OBJS = create_problem.o driver_OMP.o process_answer.o

BIN = /tmp

LU_driver : $(OBJS) $(PLAPACKLIB) 
	rm -f LU_driver*
	$(LINKER) -o $(BIN)/LU_driver.x $(LFLAGS) $(OBJS) $(PLAPACKLIB) $(LIB)
	touch LU_driver

clean :
	rm -f LU_driver* *.o core *~

