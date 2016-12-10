include ../../Make.include

OBJS = create_problem.o driver.o process_answer.o

QR_driver : $(OBJS) $(PLAPACKLIB)
	rm -f QR_driver
	$(LINKER) -o /tmp/QR_driver.x $(LFLAGS) $(OBJS) $(PLAPACKLIB) $(LIB)
	touch QR_driver

clean :
	rm -f QR_driver* *.o core *~
