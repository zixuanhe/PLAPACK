include ../Make.include

api = 	PLA_API.o  \
	PLA_API_util.o  

API:    $(api)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch API

clean:
	rm -f $(api) core *.o *~

