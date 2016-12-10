include ../../Make.include

OBJS = create_problem.o driver.o 

Chol_driver : $(OBJS) $(PLAPACKLIB) 
	rm -f Chol_driver*
	$(LINKER) -o /tmp/Chol_driver.x $(LFLAGS) $(OBJS) $(PLAPACKLIB) $(LIB)
	touch Chol_driver

clean :
	rm -f Chol_driver* *.o core *~
