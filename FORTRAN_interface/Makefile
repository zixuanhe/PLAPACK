include ../Make.include

fortran    =    PLA_Fortran.o \
	        PLA_Conversion.o \
		pla_type_conversion_f.o \
		random.o

FORTRAN_interface: $(fortran)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch FORTRAN_interface

clean:
	rm -f $(fortran) *~ core

