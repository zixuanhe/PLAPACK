include ../Make.include

TEMP = 	PLA_Temp.o \
	PLA_Environ.o \
	PLA_Init.o \
	PLA_Temp_distant_length.o

TEMP:   $(TEMP)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch TEMP

clean:
	rm -f $(TEMP) *~ core

