include ../Make.include

REDUCE = \
	PLA_Reduce.o  \
	PLA_Reduce_add.o  \
        PLA_Reduce_old.o  \
        PLA_Reduce_x.o  \
        PLA_Reduce_xx.o  \
        PLA_Perform_local_part_of_reduce_xx.o  

REDUCE: $(REDUCE)
	$(AR) vr $(PLAPACKLIB) $? 
	$(RANLIB) $(PLAPACKLIB)
	touch REDUCE

clean:
	rm -f $(REDUCE) *~ core

