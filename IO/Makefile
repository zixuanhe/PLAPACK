include ../Make.include

IO =	\
        PLA_File.o \
	PLA_File_async.o \
	PLA_File_stats.o

IO: $(IO)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch IO

clean:
	rm -f $(IO) core *.o *~

